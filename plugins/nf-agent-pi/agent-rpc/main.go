// Copyright 2013-2026, Seqera Labs
// SPDX-License-Identifier: Apache-2.0

package main

import (
	"bufio"
	"bytes"
	"context"
	"crypto/sha256"
	"crypto/tls"
	"encoding/hex"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"os/exec"
	"os/signal"
	"strings"
	"syscall"
	"time"

	"google.golang.org/grpc"
	"google.golang.org/grpc/connectivity"
	"google.golang.org/grpc/credentials"
	"google.golang.org/grpc/credentials/insecure"
	"google.golang.org/grpc/encoding"
	"google.golang.org/grpc/keepalive"
)

const connectMethod = "/nextflow.agent.AgentBroker/Connect"

type jsonCodec struct{}

func (jsonCodec) Name() string                       { return "json" }
func (jsonCodec) Marshal(v any) ([]byte, error)      { return json.Marshal(v) }
func (jsonCodec) Unmarshal(data []byte, v any) error { return json.Unmarshal(data, v) }

type options struct {
	endpoint     string
	invocationID string
	token        string
	// The decoded SHA-256 certificate digest, not the text of --fingerprint: parseArgs owns the
	// decode (see below), so nothing downstream can be handed a digest that has not been checked.
	pinned       []byte
	insecureTLS  bool
	startup      time.Duration
	connect      time.Duration
	maxLineBytes int
	harness      []string
}

type frame map[string]any

func main() {
	opts, err := parseArgs(os.Args[1:])
	if err != nil {
		fmt.Fprintln(os.Stderr, "agent-rpc:", err)
		os.Exit(2)
	}
	if err := run(opts); err != nil {
		fmt.Fprintln(os.Stderr, "agent-rpc:", err)
		os.Exit(1)
	}
}

func parseArgs(args []string) (options, error) {
	var opts options
	separator := -1
	for i, arg := range args {
		if arg == "--" {
			separator = i
			break
		}
	}
	if separator < 0 || separator == len(args)-1 {
		return opts, errors.New("expected harness command after --")
	}
	var fingerprint string
	flags := flag.NewFlagSet("agent-rpc", flag.ContinueOnError)
	flags.SetOutput(io.Discard)
	flags.StringVar(&opts.endpoint, "endpoint", "", "Nextflow driver gRPC endpoint")
	flags.StringVar(&opts.invocationID, "invocation", "", "agent invocation identity")
	flags.StringVar(&opts.token, "token", "", "single-invocation capability token")
	flags.StringVar(&fingerprint, "fingerprint", "", "SHA-256 digest of the driver TLS certificate to pin, as hex")
	flags.BoolVar(&opts.insecureTLS, "insecure", false, "dial the driver in cleartext; only valid when the driver sets agent.rpc.tls = false")
	flags.DurationVar(&opts.startup, "startup-timeout", 30*time.Second, "harness startup timeout")
	flags.DurationVar(&opts.connect, "connect-timeout", defaultConnectTimeout, "driver connection timeout")
	flags.IntVar(&opts.maxLineBytes, "max-line-bytes", 10*1024*1024, "maximum JSONL frame size")
	if err := flags.Parse(args[:separator]); err != nil {
		return opts, err
	}
	if opts.endpoint == "" || opts.invocationID == "" || opts.token == "" {
		return opts, errors.New("--endpoint, --invocation, and --token are required")
	}
	// TLS is on by default on the driver, so the pin is mandatory unless cleartext is
	// asked for explicitly. Absence of --fingerprint must never mean "dial unpinned".
	if fingerprint == "" && !opts.insecureTLS {
		return opts, errors.New("--fingerprint is required unless --insecure is given")
	}
	if fingerprint != "" && opts.insecureTLS {
		return opts, errors.New("--fingerprint and --insecure are mutually exclusive")
	}
	// The digest is decoded HERE, while we are still only reading argv, and not lazily at dial
	// time. dialCredentials is not reached until the Node harness has been forked, started, and
	// waited on for up to --startup-timeout, so a one-character typo in the digest used to cost a
	// full Pi-harness boot, report itself in .command.err underneath whatever the harness had
	// already written to stderr, and exit 1 -- the code that means "the run failed" -- rather than
	// the 2 that every other malformed argument gets. Rejecting it before anything is spawned
	// costs nothing and keeps all argument errors indistinguishable to the caller.
	//
	// Placed after the two checks above so their messages, which name the flags rather than the
	// value, stay the first thing an operator sees when the flags themselves are wrong.
	if !opts.insecureTLS {
		digest, err := parseFingerprint(fingerprint)
		if err != nil {
			return opts, err
		}
		opts.pinned = digest
	}
	if opts.maxLineBytes <= 0 {
		return opts, errors.New("--max-line-bytes must be positive")
	}
	// Zero or negative would expire the budget before the first SYN, failing every task with a
	// timeout no address could satisfy; there is no "wait forever" spelling to reserve here.
	if opts.connect <= 0 {
		return opts, errors.New("--connect-timeout must be positive")
	}
	opts.harness = append([]string(nil), args[separator+1:]...)
	return opts, nil
}

func run(opts options) error {
	ctx, cancel := signal.NotifyContext(context.Background(), os.Interrupt, syscall.SIGTERM)
	defer cancel()

	cmd := exec.CommandContext(ctx, opts.harness[0], opts.harness[1:]...)
	cmd.Stderr = os.Stderr
	childOut, err := cmd.StdoutPipe()
	if err != nil {
		return fmt.Errorf("open harness stdout: %w", err)
	}
	childIn, err := cmd.StdinPipe()
	if err != nil {
		return fmt.Errorf("open harness stdin: %w", err)
	}
	if err := cmd.Start(); err != nil {
		return fmt.Errorf("start harness: %w", err)
	}
	// Always reap the child, whatever exit path (including a panic) run() takes.
	defer terminate(cmd)

	childFrames := make(chan frame, 16)
	childErrors := make(chan error, 1)
	go scanHarness(childOut, opts.maxLineBytes, childFrames, childErrors)

	var ready frame
	select {
	case ready = <-childFrames:
		if ready["type"] != "ready" {
			return fmt.Errorf("expected harness ready frame, received %v", ready["type"])
		}
	case err := <-childErrors:
		return fmt.Errorf("harness failed before ready: %w", err)
	case <-time.After(opts.startup):
		return errors.New("timed out waiting for harness ready frame")
	case <-ctx.Done():
		return ctx.Err()
	}

	encoding.RegisterCodec(jsonCodec{})
	conn, err := grpc.NewClient(opts.endpoint, dialOptions(opts)...)
	if err != nil {
		return fmt.Errorf("connect to driver: %w", err)
	}
	defer conn.Close()
	if err := awaitDriver(ctx, conn, opts.endpoint, opts.connect); err != nil {
		return err
	}

	desc := &grpc.StreamDesc{StreamName: "Connect", ClientStreams: true, ServerStreams: true}
	stream, err := conn.NewStream(ctx, desc, connectMethod)
	if err != nil {
		return fmt.Errorf("open driver stream to %s: %w", opts.endpoint, err)
	}
	if err := stream.SendMsg(frame{"type": "connect", "invocationId": opts.invocationID, "token": opts.token}); err != nil {
		// This is the send whose refusal the operator most needs to read -- an expired or
		// already-consumed capability is answered here -- and nothing is receiving yet.
		return fmt.Errorf("register invocation: %w", sendStatus(stream, err))
	}

	hostFrames := make(chan frame, 16)
	hostErrors := make(chan error, 1)
	go recvHost(stream, hostFrames, hostErrors)
	// The ONLY writer to the harness' stdin, and it never leaves this goroutine: both Encode
	// calls below are branches of the one select loop, and the channels are what cross goroutine
	// boundaries. So no lock -- a mutex here would only claim a concurrency that does not exist.
	//
	// One encoder for the whole loop rather than one per write. That is not merely tidier: a
	// json.Encoder LATCHES its first write error and short-circuits every later Encode, which a
	// fresh per-write encoder would not. Unobservable here, because neither call outlives a
	// failure -- the hostFrames branch returns the error, and the ctx.Done branch is the last
	// write before the loop returns. Keep it that way: an Encode added below that carries on
	// after an error would silently write nothing.
	enc := json.NewEncoder(childIn)

	for {
		select {
		case msg := <-hostFrames:
			if err := enc.Encode(msg); err != nil {
				return fmt.Errorf("write harness frame: %w", err)
			}
		case msg := <-childFrames:
			msg["invocationId"] = opts.invocationID
			if err := stream.SendMsg(msg); err != nil {
				if errors.Is(err, io.EOF) {
					// Same io.EOF-hides-the-status shape as the connect send, but recvHost owns
					// RecvMsg by now, so asking for it here would be the concurrent receive grpc-go
					// forbids. It is already reading; give it a moment to hand the status over.
					select {
					case hostErr := <-hostErrors:
						return fmt.Errorf("driver stream: %w", hostErr)
					case <-time.After(2 * time.Second):
					}
				}
				return fmt.Errorf("send driver frame: %w", err)
			}
			switch msg["type"] {
			case "complete":
				if err := json.NewEncoder(os.Stdout).Encode(msg); err != nil {
					return fmt.Errorf("write final result: %w", err)
				}
				_ = stream.CloseSend()
				_ = childIn.Close()
				// Wait for the driver to END the stream before the deferred conn.Close runs.
				// CloseSend half-closes only THIS direction: the driver still has trailers to
				// flush, and dropping the TLS connection under that flush is what makes it log
				// `SSLEngine closed already` / "Transport failed" on a run that in fact
				// succeeded. recvHost owns RecvMsg, so the end of the stream reaches us as the
				// io.EOF it reports. Bounded, so a driver that never ends the stream costs this
				// much and no more.
				select {
				case <-hostErrors:
				case <-time.After(gracefulCloseTimeout):
				}
				return nil
			case "error":
				_ = stream.CloseSend()
				return fmt.Errorf("harness error [%v]: %v", msg["code"], msg["message"])
			}
		case err := <-childErrors:
			return fmt.Errorf("harness protocol: %w", err)
		case err := <-hostErrors:
			return fmt.Errorf("driver stream: %w", err)
		case <-ctx.Done():
			_ = enc.Encode(frame{"type": "cancel", "invocationId": opts.invocationID, "reason": ctx.Err().Error()})
			return ctx.Err()
		}
	}
}

// sendStatus resolves what SendMsg reported into an error that still carries the stream's status.
//
// grpc-go returns a bare io.EOF from SendMsg whenever the stream was ended by the server or the
// transport rather than by this client, and documents the status as retrievable only from RecvMsg
// ("otherwise, io.EOF is returned and the status of the stream may be discovered using RecvMsg").
// Returning that io.EOF as-is costs precisely the diagnosis the broker went to the trouble of
// sending: an UNAUTHENTICATED description reaches the task's .command.err as "EOF".
//
// Only safe to call while nothing else is receiving on the stream -- grpc-go permits concurrent
// SendMsg and RecvMsg, but not two concurrent RecvMsg. Once recvHost is running it owns the
// receive side, so the frame loop recovers its status from hostErrors instead.
func sendStatus(stream grpc.ClientStream, err error) error {
	if !errors.Is(err, io.EOF) {
		// Anything else was generated by this client and is already the status it wants reported.
		return err
	}
	var discard frame
	if recvErr := stream.RecvMsg(&discard); recvErr != nil && !errors.Is(recvErr, io.EOF) {
		return recvErr
	}
	// A clean close on both sides leaves nothing better to say, and callers still test for io.EOF.
	return err
}

// How long the driver connection has to come up before the task is failed.
//
// The driver no longer always hands this proxy an address a human wrote: where the endpoint used to
// be either a container engine's host alias or an explicit agent.rpc.remoteHost, it is now usually
// INFERRED from the driver's own default route, and the expensive way for inference to be wrong is
// to be plausible -- an address that routes nowhere rather than one that answers "refused".
//
// 30s, matching --startup-timeout, is above grpc-go's own 20s minimum per-attempt connect timeout,
// so a slow TCP+TLS handshake over a loaded link still completes inside the budget and only an
// endpoint that never answers spends it.
const defaultConnectTimeout = 30 * time.Second

// awaitDriver bounds the wait for the driver connection to become usable, so a wrong endpoint fails
// the task in seconds and names the address it tried.
//
// A REFUSED endpoint is already loud without this: grpc.NewClient is lazy, the stream below is
// opened without WaitForReady, and a channel in TRANSIENT_FAILURE fails NewStream immediately. A
// BLACKHOLED one -- the shape a mis-inferred address takes, packets silently dropped -- is not:
// there is nothing to fail on, the channel sits in CONNECTING, and both NewStream and the OS would
// wait on it. Meanwhile the driver holds this invocation's capability until agent.rpc.capabilityTimeout
// expires, an hour by default, pinning the request behind a connection that will never arrive.
//
// TRANSIENT_FAILURE returns nil rather than an error of its own, deliberately: this state carries no
// reason, while the status NewStream fails with quotes the transport's ("connection refused", or the
// pinning refusal from pinnedTLSConfig). Reporting it here would cost exactly the diagnosis the
// operator needs. grpc's initial reconnect backoff is a second, so the caller is in NewStream long
// before the channel leaves that state.
//
// Nothing is offered to the endpoint here. The capability token crosses the wire in run()'s connect
// frame, after this returns and after the stream opened, which is after the TLS pin has verified --
// so a driver that fails the pin is never told the token.
func awaitDriver(ctx context.Context, conn *grpc.ClientConn, endpoint string, budget time.Duration) error {
	dialCtx, cancel := context.WithTimeout(ctx, budget)
	defer cancel()
	// NewClient connects nothing until a stream asks for a subchannel, so without this the wait
	// below would time out against a channel that never left IDLE.
	conn.Connect()
	for {
		state := conn.GetState()
		if state == connectivity.Ready || state == connectivity.TransientFailure {
			return nil
		}
		if !conn.WaitForStateChange(dialCtx, state) {
			if err := ctx.Err(); err != nil {
				// The run was cancelled or the task killed; the endpoint is not at fault.
				return err
			}
			return fmt.Errorf("connect to driver at %s: no connection after %s", endpoint, budget)
		}
	}
}

// parseFingerprint decodes the pinned SHA-256 certificate digest. The driver emits lowercase
// unseparated hex over the certificate's DER; colons and case are tolerated here only so a
// digest copied out of keytool or `openssl x509 -fingerprint` also works.
func parseFingerprint(value string) ([]byte, error) {
	cleaned := strings.ReplaceAll(strings.ToLower(strings.TrimSpace(value)), ":", "")
	digest, err := hex.DecodeString(cleaned)
	if err != nil {
		return nil, fmt.Errorf("invalid --fingerprint %q: expected a hex SHA-256 digest", value)
	}
	if len(digest) != sha256.Size {
		return nil, fmt.Errorf("invalid --fingerprint %q: expected %d hex bytes, got %d", value, sha256.Size, len(digest))
	}
	return digest, nil
}

// How often the proxy pings an otherwise silent driver connection, and how long the driver then
// has to answer. The post-connect capability deadline is deliberately gone -- an agent may run
// for hours -- which leaves nothing on THIS side watching a driver that vanished: a killed JVM
// or a lost driver node sends no FIN, so the blocking RecvMsg in recvHost would simply never
// return and the task would burn its container until the operating system's own TCP keepalive
// noticed, 2h on Linux.
//
// Time is 10 minutes and is NOT derived from the 10s the current broker permits. This binary
// ships in the runner container image (/usr/local/bin/agent-rpc, see PiAgentRunner), not in the
// plugin jar, so it routinely meets a driver it was not built alongside -- including one whose
// broker never calls permitKeepAliveTime and therefore enforces grpc-java's default of 5 minutes.
// Two pings inside that window is a GOAWAY, and a GOAWAY is unrecoverable here: the capability
// token is single-use and was consumed on connect, so grpc-go's reconnect cannot register again
// and the invocation is lost. 10 minutes leaves 2x margin on the tightest policy any driver
// enforces; 5 would be arithmetically equal to it, so ordinary clock skew alone would turn a
// legitimate ping into a strike. The cost of erring long is only detection latency, and even that
// is mostly theoretical: the driver pings every 60s and grpc-go re-arms this timer from the last
// READ, so against a live driver this ping is essentially never sent.
//
// Timeout is 90 seconds and does double duty, which is the part that is easy to miss. Setting
// Time at all makes grpc-go call setsockopt(TCP_USER_TIMEOUT, Timeout) on the socket -- a real
// kernel timer on Linux, where these tasks run -- so this value is also the budget for
// unacknowledged OUTBOUND data, replacing the ~15 minutes the default tcp_retries2 allows. It is
// set well above grpc-go's own 20s default so a network blip or a stop-the-world GC on the driver
// is not fatal, and deliberately just above the ~80s (keepAliveTime 60s + keepAliveTimeout 20s,
// see AgentRpcBroker) at which the driver's server keepalive already gives up on this task: past
// that point the driver has failed the task anyway, so nothing is lost by the client's kernel
// giving up too, and below it the kernel could kill a link the driver still considers good.
//
// PermitWithoutStream stays false, matching the broker's permitKeepAliveWithoutCalls. grpc-go
// then parks the keepalive goroutine while no stream is open, which is exactly right: a ping sent
// before the Connect stream exists would be a strike against a server that is not obliged to
// tolerate it, and the window it would cover is the sub-millisecond gap between dial and stream.
var keepaliveParams = keepalive.ClientParameters{
	Time:                10 * time.Minute,
	Timeout:             90 * time.Second,
	PermitWithoutStream: false,
}

// dialOptions is the whole transport configuration for the driver connection, kept in one place
// and shared with the tests so they exercise the transport the binary actually dials with rather
// than a hand-assembled approximation of it.
func dialOptions(opts options) []grpc.DialOption {
	return []grpc.DialOption{
		grpc.WithTransportCredentials(dialCredentials(opts)),
		grpc.WithDefaultCallOptions(grpc.ForceCodec(jsonCodec{})),
		grpc.WithKeepaliveParams(keepaliveParams),
	}
}

// dialCredentials cannot fail: parseArgs has already decoded and length-checked the pin, so the
// only two outcomes are the deliberate cleartext escape hatch and a pinned TLS transport.
func dialCredentials(opts options) credentials.TransportCredentials {
	if opts.insecureTLS {
		return insecure.NewCredentials()
	}
	return credentials.NewTLS(pinnedTLSConfig(opts.pinned))
}

// pinnedTLSConfig pins the driver's per-run, self-signed certificate by SHA-256 digest. There is
// no CA to chain to and no name to match, so the digest is the whole trust decision:
// InsecureSkipVerify turns off chain and hostname verification only, and the hook below still
// runs on every handshake.
//
// The digest is taken over the leaf certificate's DER encoding, which is what the driver hashes
// too (X509Certificate.getEncoded(), see AgentRpcTlsCredentials). Hashing the PEM text or its
// base64 body instead yields an equally well-formed digest that simply never matches, and the
// failure is indistinguishable from a genuine pinning rejection.
//
// The comparison lives in VerifyConnection and NOT in VerifyPeerCertificate, and moving it back
// silently unpins any resumed connection. crypto/tls does not re-verify certificates on a
// resumption: handshake_client_tls13.go returns from readServerCertificate the moment hs.usingPSK
// is set, and the TLS 1.2 session-ticket path never reaches verifyServerCertificate at all -- so
// VerifyPeerCertificate is not consulted on either, while VerifyConnection is called on both.
// That cannot fire today (grpc-go leaves ClientSessionCache nil, so no PSK is ever offered, and
// each agent-rpc process opens exactly one connection), which is precisely the danger: a session
// cache or a shared channel introduced later would turn the pin off with nothing to notice. It is
// the same number of lines either way, so there is no reason to sit on the fragile one.
//
// cs.PeerCertificates[0].Raw carries the same DER that rawCerts[0] would have: on a full
// handshake it is the input crypto/tls handed to x509.ParseCertificate, and on a resumed one it
// is restored from the session that the full handshake already pinned.
//
// A nil pin fails closed rather than open -- bytes.Equal cannot match a 32-byte digest against
// nothing -- so a driver-side regression that stopped emitting a fingerprint would refuse every
// certificate instead of accepting every certificate.
func pinnedTLSConfig(pinned []byte) *tls.Config {
	return &tls.Config{
		InsecureSkipVerify: true,
		VerifyConnection: func(cs tls.ConnectionState) error {
			if len(cs.PeerCertificates) == 0 {
				return errors.New("driver presented no TLS certificate")
			}
			presented := sha256.Sum256(cs.PeerCertificates[0].Raw)
			if !bytes.Equal(presented[:], pinned) {
				// A fingerprint is a public commitment, not a secret, so both values are
				// safe to report and make the failure unambiguously a pinning failure.
				return fmt.Errorf("driver TLS certificate fingerprint mismatch: pinned %s, presented %s",
					hex.EncodeToString(pinned), hex.EncodeToString(presented[:]))
			}
			return nil
		},
	}
}

func scanHarness(source io.Reader, max int, output chan<- frame, failures chan<- error) {
	defer func() {
		if r := recover(); r != nil {
			failures <- fmt.Errorf("harness scan panic: %v", r)
		}
	}()
	scanner := bufio.NewScanner(source)
	scanner.Buffer(make([]byte, 64*1024), max)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		var msg frame
		if err := json.Unmarshal([]byte(line), &msg); err != nil {
			failures <- fmt.Errorf("malformed JSON: %w", err)
			return
		}
		// A literal `null` (or any JSON null) unmarshals into a nil map with no
		// error; forwarding it would panic when the bridge writes into the map.
		if msg == nil {
			continue
		}
		output <- msg
	}
	if err := scanner.Err(); err != nil {
		failures <- err
		return
	}
	failures <- io.EOF
}

// How long to wait, after the final frame, for the driver to end the stream. Only the tail of a
// successful invocation waits: the result has already been written to stdout by then, so this costs
// nothing but a clean connection teardown.
const gracefulCloseTimeout = 2 * time.Second

func recvHost(stream grpc.ClientStream, output chan<- frame, failures chan<- error) {
	defer func() {
		if r := recover(); r != nil {
			failures <- fmt.Errorf("host receive panic: %v", r)
		}
	}()
	for {
		var msg frame
		if err := stream.RecvMsg(&msg); err != nil {
			failures <- err
			return
		}
		output <- msg
	}
}

func terminate(cmd *exec.Cmd) {
	if cmd == nil || cmd.Process == nil {
		return
	}
	_ = cmd.Process.Signal(syscall.SIGTERM)
	done := make(chan struct{})
	go func() { _ = cmd.Wait(); close(done) }()
	select {
	case <-done:
	case <-time.After(2 * time.Second):
		_ = cmd.Process.Kill()
		<-done
	}
}
