// Copyright 2013-2026, Seqera Labs
// SPDX-License-Identifier: Apache-2.0

package main

import (
	"context"
	"crypto/ecdsa"
	"crypto/elliptic"
	"crypto/rand"
	"crypto/sha256"
	"crypto/tls"
	"crypto/x509"
	"crypto/x509/pkix"
	"encoding/hex"
	"errors"
	"io"
	"math/big"
	"net"
	"strings"
	"testing"
	"time"

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/credentials"
	"google.golang.org/grpc/status"
)

func TestParseArgs(t *testing.T) {
	opts, err := parseArgs([]string{
		"--endpoint", "driver:1234",
		"--invocation", "inv-1",
		"--token", "secret",
		"--fingerprint", strings.Repeat("ab", sha256.Size),
		"--startup-timeout", "5s",
		"--", "node", "runner.mjs",
	})
	if err != nil {
		t.Fatal(err)
	}
	if opts.endpoint != "driver:1234" || opts.invocationID != "inv-1" || opts.token != "secret" {
		t.Fatalf("unexpected options: %#v", opts)
	}
	if hex.EncodeToString(opts.pinned) != strings.Repeat("ab", sha256.Size) || opts.insecureTLS {
		t.Fatalf("unexpected transport options: %#v", opts)
	}
	if opts.startup != 5*time.Second || strings.Join(opts.harness, " ") != "node runner.mjs" {
		t.Fatalf("unexpected runtime options: %#v", opts)
	}
	// Not passed above: the driver does not emit --connect-timeout, so the default is what every
	// task actually runs with.
	if opts.connect != defaultConnectTimeout {
		t.Fatalf("expected the shipped connect budget by default, got %v", opts.connect)
	}
}

func TestParseArgsRejectsANonPositiveConnectTimeout(t *testing.T) {
	args := []string{
		"--endpoint", "driver:1234", "--invocation", "inv-1", "--token", "secret",
		"--insecure", "--connect-timeout", "0", "--", "node",
	}
	if _, err := parseArgs(args); err == nil {
		t.Fatal("expected a zero connect budget to be rejected: it would fail every task before the first SYN")
	}
}

func TestParseArgsRejectsMissingIdentity(t *testing.T) {
	if _, err := parseArgs([]string{"--endpoint", "driver:1234", "--insecure", "--", "node"}); err == nil {
		t.Fatal("expected missing invocation identity to fail")
	}
}

// The driver serves TLS unless explicitly told not to, so a missing pin must fail closed: if
// absence meant "dial unpinned", a driver-side regression would silently downgrade the transport.
func TestParseArgsRequiresFingerprintUnlessInsecure(t *testing.T) {
	base := []string{"--endpoint", "driver:1234", "--invocation", "inv-1", "--token", "secret"}
	if _, err := parseArgs(append(append([]string(nil), base...), "--", "node")); err == nil {
		t.Fatal("expected a missing --fingerprint to fail")
	}
	if _, err := parseArgs(append(append([]string(nil), base...), "--insecure", "--", "node")); err != nil {
		t.Fatalf("expected --insecure to stand in for the pin: %v", err)
	}
}

func TestParseArgsRejectsFingerprintWithInsecure(t *testing.T) {
	args := []string{
		"--endpoint", "driver:1234", "--invocation", "inv-1", "--token", "secret",
		"--fingerprint", strings.Repeat("ab", sha256.Size), "--insecure", "--", "node",
	}
	if _, err := parseArgs(args); err == nil {
		t.Fatal("expected --fingerprint and --insecure to be mutually exclusive")
	}
}

func TestParseFingerprintForms(t *testing.T) {
	canonical := strings.Repeat("ab", sha256.Size)
	// Lowercase unseparated hex is what the driver emits; the rest are tolerated so a digest
	// pasted from keytool or `openssl x509 -fingerprint` also works.
	for _, in := range []string{canonical, strings.ToUpper(canonical), colonize(strings.ToUpper(canonical)), "  " + canonical + "  "} {
		got, err := parseFingerprint(in)
		if err != nil {
			t.Fatalf("%q: %v", in, err)
		}
		if hex.EncodeToString(got) != canonical {
			t.Fatalf("%q decoded to %s", in, hex.EncodeToString(got))
		}
	}
	// "sha256:" prefixes only ever precede base64 in the SSH form, never hex, so accepting one
	// here would accept a string nothing produces.
	for _, bad := range []string{"", "zz", "not-hex", canonical + "ab", canonical[:len(canonical)-2], "sha256:" + canonical} {
		if _, err := parseFingerprint(bad); err == nil {
			t.Fatalf("expected %q to be rejected", bad)
		}
	}
}

// A malformed digest has to be caught while argv is still the only thing that has happened. By
// the time dialCredentials runs, the Node harness has been forked and waited on for up to
// --startup-timeout, so a typo used to cost a full harness boot and then report itself as a run
// failure (exit 1) rather than as the bad argument it is (exit 2).
func TestParseArgsRejectsMalformedFingerprint(t *testing.T) {
	args := []string{
		"--endpoint", "driver:1234", "--invocation", "inv-1", "--token", "secret",
		"--fingerprint", "not-a-digest", "--", "node",
	}
	_, err := parseArgs(args)
	if err == nil {
		t.Fatal("expected a malformed pin to be rejected while parsing arguments")
	}
	if !strings.Contains(err.Error(), "--fingerprint") {
		t.Fatalf("the error should name the flag at fault: %v", err)
	}
}

func TestDialCredentialsInsecureOptsOutOfTLS(t *testing.T) {
	if got := dialCredentials(options{insecureTLS: true}).Info().SecurityProtocol; got != "insecure" {
		t.Fatalf("expected the escape hatch to select cleartext, got %q", got)
	}
}

func TestPinnedFingerprintAcceptsTheDriverCertificate(t *testing.T) {
	cert, fingerprint := selfSignedCertificate(t)
	if err := dialPinned(t, serveTLS(t, cert), fingerprint); err != nil {
		t.Fatalf("expected the handshake to succeed with a matching pin: %v", err)
	}
}

// dialCredentials can no longer report an error, so the fail-closed property rests entirely on
// bytes.Equal refusing to match a 32-byte digest against nothing. If that were ever softened into
// "no pin means accept anything", a driver-side regression that stopped emitting --fingerprint
// would unpin every task in the run without a single failure to show for it.
func TestDialCredentialsWithoutAPinRefusesEveryCertificate(t *testing.T) {
	cert, _ := selfSignedCertificate(t)
	err := dialThrough(t, serveTLS(t, cert), options{})
	if err == nil {
		t.Fatal("expected an absent pin to refuse the certificate, not to accept it")
	}
	if !strings.Contains(err.Error(), "fingerprint mismatch") {
		t.Fatalf("expected a pinning refusal, got: %v", err)
	}
}

// TLS session resumption is where a certificate pin most often goes quietly missing: crypto/tls
// does not re-verify certificates on a resumed handshake, so a check installed in
// VerifyPeerCertificate is simply never consulted. grpc-go opens no session cache today, which is
// why this has to reach the TLS config directly rather than dial through grpc -- nothing in the
// shipped path can currently produce a resumption, and therefore nothing in the shipped path
// would notice the day a cache or a shared channel appears. What is asserted here is the
// placement of our own hook, which is the entire change.
func TestPinSurvivesASessionResumption(t *testing.T) {
	cert, served := selfSignedCertificate(t)
	endpoint := serveRawTLS(t, cert)
	cache := tls.NewLRUClientSessionCache(4)

	resumed, err := pinnedHandshake(t, endpoint, served, cache)
	if err != nil {
		t.Fatalf("first handshake: %v", err)
	}
	if resumed {
		t.Fatal("the first handshake cannot be a resumption")
	}
	// Proven, not assumed: if the second handshake were full, the rest of this test would pass
	// against the very placement it exists to reject.
	if resumed, err = pinnedHandshake(t, endpoint, served, cache); err != nil {
		t.Fatalf("second handshake: %v", err)
	} else if !resumed {
		t.Fatal("expected the second handshake to resume the session")
	}

	// Same warm cache, so the same PSK is offered and no certificate crosses the wire -- but the
	// identity behind it is still the server's, and it still has to match.
	_, other := selfSignedCertificate(t)
	if _, err := pinnedHandshake(t, endpoint, other, cache); err == nil {
		t.Fatal("a resumed handshake was accepted without checking the pin")
	} else if !strings.Contains(err.Error(), "fingerprint mismatch") {
		t.Fatalf("expected a pinning refusal on the resumed handshake, got: %v", err)
	}
}

// The values, and the fact that they reach the dial at all. An end-to-end ping exchange at the
// shipped interval is not testable in a unit suite (and testing it at an interval we do not ship
// would assert nothing about what we do ship), so this pins the decision instead.
func TestKeepaliveStaysWithinTheOldestDriversPingBudget(t *testing.T) {
	// grpc-java's ServerBuilder default when a broker never calls permitKeepAliveTime. This
	// binary lives in the runner container image, not in the plugin jar, so it meets such a
	// driver whenever an image and a driver are not upgraded together.
	const grpcJavaDefaultPermitKeepAliveTime = 5 * time.Minute
	if keepaliveParams.Time < 2*grpcJavaDefaultPermitKeepAliveTime {
		// Not merely "greater than": at exactly the enforcement floor, ordinary clock skew turns
		// a legitimate ping into a strike, two strikes are a GOAWAY, and a GOAWAY cannot be
		// recovered from because the capability token was consumed on connect.
		t.Fatalf("ping interval %v leaves no margin over the oldest driver's %v enforcement floor",
			keepaliveParams.Time, grpcJavaDefaultPermitKeepAliveTime)
	}
	// Timeout is also handed to setsockopt(TCP_USER_TIMEOUT) by grpc-go, so it is simultaneously
	// the kernel's budget for unacknowledged outbound data. Below the ~80s at which the driver's
	// own server keepalive gives up (AgentRpcBroker: 60s + 20s), it would start killing links the
	// driver still considers healthy -- a far tighter policy than anything decided here.
	if keepaliveParams.Timeout <= 80*time.Second {
		t.Fatalf("ack budget %v undercuts the driver's own ~80s liveness window", keepaliveParams.Timeout)
	}
	// Must match the broker's permitKeepAliveWithoutCalls, which is deliberately left false.
	if keepaliveParams.PermitWithoutStream {
		t.Fatal("pinging with no stream open is a strike against a driver that does not permit it")
	}
	// Cheap tripwire for the one thing no behavioural test can observe: that the parameters are
	// actually on the dial. Credentials, codec, keepalive.
	if got := len(dialOptions(options{pinned: make([]byte, sha256.Size)})); got != 3 {
		t.Fatalf("expected credentials, codec and keepalive on the dial, got %d options", got)
	}
}

// The digest, not the certificate's issuer or names, is the whole trust decision, so a server
// this proxy did not expect must be refused with an unmistakable pinning error rather than a
// generic handshake failure an operator cannot act on.
func TestPinnedFingerprintRejectsAnotherCertificate(t *testing.T) {
	cert, served := selfSignedCertificate(t)
	_, other := selfSignedCertificate(t)
	err := dialPinned(t, serveTLS(t, cert), other)
	if err == nil {
		t.Fatal("expected a certificate that does not match the pin to be refused")
	}
	if !strings.Contains(err.Error(), "fingerprint mismatch") {
		t.Fatalf("pinning error lost on the way out of grpc: %v", err)
	}
	if !strings.Contains(err.Error(), served) || !strings.Contains(err.Error(), other) {
		t.Fatalf("both the pinned and the presented digest should be reported: %v", err)
	}
}

// The endpoint is no longer always an address a human wrote: it is usually inferred from the
// driver's default route, and the expensive way for inference to be wrong is to be plausible --
// an address that routes nowhere rather than one that answers "refused". Left to itself the channel
// sits in CONNECTING, so nothing here fails, while the driver holds this invocation's capability
// for the full agent.rpc.capabilityTimeout (an hour by default) waiting for a connection that will
// never arrive.
//
// The dialer neither answers nor refuses, which is what a dropped SYN looks like from this side and
// is the one case no real endpoint can be made to reproduce on demand.
func TestAwaitDriverGivesUpOnABlackholedEndpoint(t *testing.T) {
	const endpoint = "127.0.0.1:1"
	released := make(chan struct{})
	t.Cleanup(func() { close(released) })
	blackhole := grpc.WithContextDialer(func(ctx context.Context, _ string) (net.Conn, error) {
		select {
		case <-ctx.Done():
		case <-released:
		}
		return nil, errors.New("blackholed")
	})
	conn, err := grpc.NewClient(endpoint, append(dialOptions(options{insecureTLS: true}), blackhole)...)
	if err != nil {
		t.Fatal(err)
	}
	defer conn.Close()

	started := time.Now()
	err = awaitDriver(context.Background(), conn, endpoint, 200*time.Millisecond)

	if err == nil {
		t.Fatal("expected a blackholed endpoint to fail the task, not to be waited on")
	}
	// Bounded by OUR budget, not by grpc-go's 20s per-attempt connect timeout or the kernel's own
	// retransmission schedule, either of which would dominate if the deadline were not enforced here.
	if elapsed := time.Since(started); elapsed > 5*time.Second {
		t.Fatalf("gave up only after %v, so the budget is not what bounded the wait", elapsed)
	}
	if !strings.Contains(err.Error(), endpoint) {
		t.Fatalf("the failure must name the address it tried: %v", err)
	}
}

// A refused endpoint was already fast, and has to stay that way -- but through the whole of run(),
// because naming the address is half the point: an operator reading .command.err has to be able to
// tell WHICH driver address the task tried when the inference picked the wrong interface.
func TestRunFailsFastOnARefusedEndpointAndNamesIt(t *testing.T) {
	endpoint := closedPort(t)
	opts := options{
		endpoint:     endpoint,
		invocationID: "inv-1",
		token:        "secret",
		insecureTLS:  true,
		startup:      10 * time.Second,
		connect:      10 * time.Second,
		maxLineBytes: 64 * 1024,
		// The proxy dials only after the harness is up, so a stand-in has to announce itself and
		// then stay alive; `cat` blocks on stdin exactly as the real harness does between frames.
		harness: []string{"sh", "-c", `printf '{"type":"ready"}\n'; exec cat`},
	}

	started := time.Now()
	err := run(opts)

	if err == nil {
		t.Fatal("expected a refused endpoint to fail the task")
	}
	if elapsed := time.Since(started); elapsed > opts.connect {
		t.Fatalf("a refusal took %v: it must not be waited out like a blackhole", elapsed)
	}
	if !strings.Contains(err.Error(), endpoint) {
		t.Fatalf("the failure must name the address it tried: %v", err)
	}
	// The transport's own reason survives the connect wait; awaitDriver deliberately reports
	// nothing of its own for TRANSIENT_FAILURE, which carries no reason to report.
	if !strings.Contains(err.Error(), "refused") {
		t.Fatalf("the transport's reason was lost on the way out: %v", err)
	}
}

// closedPort is an address nothing listens on: bound and then released, so the port is a real free
// one rather than a guess that something else on the machine might occupy.
func closedPort(t *testing.T) string {
	t.Helper()
	listener, err := net.Listen("tcp", "127.0.0.1:0")
	if err != nil {
		t.Fatal(err)
	}
	endpoint := listener.Addr().String()
	if err := listener.Close(); err != nil {
		t.Fatal(err)
	}
	return endpoint
}

// selfSignedCertificate stands in for the driver's per-run identity, and returns the fingerprint
// in the driver's own form: SHA-256 over the certificate DER, lowercase hex, no separators.
//
// It mirrors what AgentRpcTlsCredentials actually emits -- EC P-256, self-signed, the same two
// names, and deliberately no basicConstraints and no key-usage extensions -- so the pin path is
// exercised against the certificate shape the broker really serves rather than a tidier one.
func selfSignedCertificate(t *testing.T) (tls.Certificate, string) {
	t.Helper()
	key, err := ecdsa.GenerateKey(elliptic.P256(), rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	template := &x509.Certificate{
		SerialNumber: big.NewInt(time.Now().UnixNano()),
		Subject:      pkix.Name{CommonName: "nextflow-agent-rpc"},
		NotBefore:    time.Now().Add(-time.Hour),
		NotAfter:     time.Now().Add(24 * time.Hour),
		DNSNames:     []string{"localhost"},
		IPAddresses:  []net.IP{net.ParseIP("127.0.0.1")},
	}
	der, err := x509.CreateCertificate(rand.Reader, template, template, &key.PublicKey, key)
	if err != nil {
		t.Fatal(err)
	}
	digest := sha256.Sum256(der)
	return tls.Certificate{Certificate: [][]byte{der}, PrivateKey: key}, hex.EncodeToString(digest[:])
}

func serveTLS(t *testing.T, cert tls.Certificate) string {
	t.Helper()
	listener, err := net.Listen("tcp", "127.0.0.1:0")
	if err != nil {
		t.Fatal(err)
	}
	server := grpc.NewServer(grpc.Creds(credentials.NewTLS(&tls.Config{Certificates: []tls.Certificate{cert}})))
	go func() { _ = server.Serve(listener) }()
	t.Cleanup(server.Stop)
	return listener.Addr().String()
}

// serveRawTLS is the resumption test's server: crypto/tls directly rather than grpc, because
// grpc-go offers no way to hand the server a session-ticket policy and the whole point is to
// obtain a resumable session. It writes a byte and then blocks, so a client read has something
// to return once the post-handshake NewSessionTicket has been processed.
func serveRawTLS(t *testing.T, cert tls.Certificate) string {
	t.Helper()
	listener, err := tls.Listen("tcp", "127.0.0.1:0", &tls.Config{
		Certificates: []tls.Certificate{cert},
		MinVersion:   tls.VersionTLS13,
	})
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { _ = listener.Close() })
	go func() {
		for {
			conn, err := listener.Accept()
			if err != nil {
				return
			}
			go func() {
				defer conn.Close()
				if _, err := conn.Write([]byte{'x'}); err != nil {
					return
				}
				_, _ = io.Copy(io.Discard, conn)
			}()
		}
	}()
	return listener.Addr().String()
}

// pinnedHandshake runs one real handshake using the production pinning config with nothing added
// but the shared session cache, and reports whether it resumed. The read is not incidental: in
// TLS 1.3 the NewSessionTicket arrives after the handshake completes and the client only
// processes it while reading, so without it there would be nothing in the cache to resume from.
func pinnedHandshake(t *testing.T, endpoint, fingerprint string, cache tls.ClientSessionCache) (bool, error) {
	t.Helper()
	config := pinnedTLSConfig(mustPin(t, fingerprint))
	config.ClientSessionCache = cache
	config.MinVersion = tls.VersionTLS13
	conn, err := tls.Dial("tcp", endpoint, config)
	if err != nil {
		return false, err
	}
	defer conn.Close()
	if _, err := conn.Read(make([]byte, 1)); err != nil {
		return false, err
	}
	return conn.ConnectionState().DidResume, nil
}

// dialPinned exercises the real transport: grpc.NewClient is lazy, so the handshake -- and hence
// the pin check -- only happens once a stream is opened.
func dialPinned(t *testing.T, endpoint, fingerprint string) error {
	t.Helper()
	return dialThrough(t, endpoint, options{pinned: mustPin(t, fingerprint)})
}

// dialThrough builds the connection from dialOptions, i.e. from the exact option set run() uses,
// so these tests cannot drift away from what the binary dials with. That also means they inherit
// the JSON codec and the keepalive parameters; harmless here, since no message is ever sent and
// no test runs for anything close to the ping interval.
//
// It waits through awaitDriver for the same reason: that is the sequence run() performs, and a
// connect wait that swallowed the transport's own failure reason would show up here as a pinning
// test that stopped seeing pinning errors.
func dialThrough(t *testing.T, endpoint string, opts options) error {
	t.Helper()
	conn, err := grpc.NewClient(endpoint, dialOptions(opts)...)
	if err != nil {
		return err
	}
	defer conn.Close()
	ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
	defer cancel()
	if err := awaitDriver(ctx, conn, endpoint, 10*time.Second); err != nil {
		return err
	}
	desc := &grpc.StreamDesc{StreamName: "Connect", ClientStreams: true, ServerStreams: true}
	_, err = conn.NewStream(ctx, desc, connectMethod)
	return err
}

func mustPin(t *testing.T, fingerprint string) []byte {
	t.Helper()
	digest, err := parseFingerprint(fingerprint)
	if err != nil {
		t.Fatal(err)
	}
	return digest
}

func colonize(s string) string {
	var pairs []string
	for i := 0; i+1 < len(s); i += 2 {
		pairs = append(pairs, s[i:i+2])
	}
	return strings.Join(pairs, ":")
}

func TestScanHarness(t *testing.T) {
	frames := make(chan frame, 2)
	failures := make(chan error, 1)
	go scanHarness(strings.NewReader("\n{\"type\":\"ready\"}\n"), 1024, frames, failures)

	if got := <-frames; got["type"] != "ready" {
		t.Fatalf("unexpected frame: %#v", got)
	}
	if err := <-failures; err == nil || err.Error() != "EOF" {
		t.Fatalf("expected EOF, got %v", err)
	}
}

func TestScanHarnessRejectsMalformedJSON(t *testing.T) {
	frames := make(chan frame, 1)
	failures := make(chan error, 1)
	go scanHarness(strings.NewReader("not-json\n"), 1024, frames, failures)
	if err := <-failures; err == nil || !strings.Contains(err.Error(), "malformed JSON") {
		t.Fatalf("unexpected error: %v", err)
	}
}

// A stream whose RecvMsg answers with a status, standing in for a transport-terminated send.
//
// grpc-go returns a bare io.EOF from SendMsg whenever the stream was ended by the server or the
// transport rather than by this client, and documents the status as retrievable only from RecvMsg.
// Provoking that against a live server means killing the transport inside the window between
// NewStream and SendMsg, which is inherently racy. The contract under test -- "on io.EOF, ask
// RecvMsg" -- is exact, so it is asserted against grpc's own ClientStream interface rather than
// against a coin flip. The embedded nil interface supplies the methods this path never calls.
type statusStream struct {
	grpc.ClientStream
	recvErr   error
	recvCalls int
}

func (s *statusStream) RecvMsg(any) error {
	s.recvCalls++
	return s.recvErr
}

func TestSendStatusRecoversTheStatusBehindAnEOF(t *testing.T) {
	refusal := status.Error(codes.Unauthenticated, "Agent RPC invocation capability expired after 3600s while the task waited to start")
	stream := &statusStream{recvErr: refusal}

	err := sendStatus(stream, io.EOF)

	if stream.recvCalls != 1 {
		t.Fatalf("expected exactly one RecvMsg to recover the status, got %d", stream.recvCalls)
	}
	if status.Code(err) != codes.Unauthenticated {
		t.Fatalf("status code lost: got %v, want %v", status.Code(err), codes.Unauthenticated)
	}
	// The description is the whole point: without it the task's .command.err reads "EOF" and the
	// operator cannot tell a lapsed capability from a forged one.
	if !strings.Contains(err.Error(), "capability expired after 3600s") {
		t.Fatalf("status description lost: %v", err)
	}
}

func TestSendStatusLeavesAClientGeneratedErrorAlone(t *testing.T) {
	direct := status.Error(codes.ResourceExhausted, "grpc: trying to send message larger than max")
	stream := &statusStream{recvErr: status.Error(codes.Unavailable, "should not be consulted")}

	err := sendStatus(stream, direct)

	if err != direct {
		t.Fatalf("a client-generated status must pass through untouched: got %v", err)
	}
	if stream.recvCalls != 0 {
		t.Fatalf("RecvMsg must not be consulted when SendMsg already carried the status, got %d calls", stream.recvCalls)
	}
}

func TestSendStatusKeepsEOFWhenTheStreamCarriesNoStatus(t *testing.T) {
	stream := &statusStream{recvErr: io.EOF}

	err := sendStatus(stream, io.EOF)

	// A clean half-close on both sides leaves nothing better to report, and callers up the stack
	// still test this with errors.Is(err, io.EOF).
	if !errors.Is(err, io.EOF) {
		t.Fatalf("expected the original io.EOF to survive, got %v", err)
	}
}
