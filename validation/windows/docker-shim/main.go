/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// Command docker is a stand-in for the Docker Desktop CLI, used to test Nextflow
// on GitHub-hosted Windows runners, which cannot run Docker Desktop.
//
// It forwards every command to a Docker Engine running in a WSL2 distribution
// (the NXF_WSL_DISTRO environment variable, default "nxf"), translating the
// Windows host paths of bind mounts to their WSL location (C:\dir -> /mnt/c/dir),
// which is what Docker Desktop does for Linux containers.
package main

import (
	"errors"
	"fmt"
	"os"
	"os/exec"
	"regexp"
	"strings"
)

var windowsMount = regexp.MustCompile(`^([A-Za-z]):([\\/][^:]*)?(:.*)?$`)

// translateVolume maps the host part of a bind mount spec, e.g.
// `C:\work\ab:/c/work/ab:ro` -> `/mnt/c/work/ab:/c/work/ab:ro`
func translateVolume(spec string) string {
	m := windowsMount.FindStringSubmatch(spec)
	if m == nil {
		return spec
	}
	host := "/mnt/" + strings.ToLower(m[1]) + strings.ReplaceAll(m[2], `\`, "/")
	return host + m[3]
}

// translateArgs translates the bind mounts in the docker command line arguments
func translateArgs(args []string) []string {
	result := make([]string, len(args))
	copy(result, args)
	for i := 0; i < len(result); i++ {
		arg := result[i]
		switch {
		case (arg == "-v" || arg == "--volume") && i+1 < len(result):
			result[i+1] = translateVolume(result[i+1])
			i++
		case strings.HasPrefix(arg, "--volume="):
			result[i] = "--volume=" + translateVolume(strings.TrimPrefix(arg, "--volume="))
		}
	}
	return result
}

func main() {
	distro := os.Getenv("NXF_WSL_DISTRO")
	if distro == "" {
		distro = "nxf"
	}
	args := append([]string{"-d", distro, "-u", "root", "--exec", "docker"}, translateArgs(os.Args[1:])...)
	cmd := exec.Command("wsl.exe", args...)
	cmd.Stdin = os.Stdin
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		var exit *exec.ExitError
		if errors.As(err, &exit) {
			os.Exit(exit.ExitCode())
		}
		fmt.Fprintln(os.Stderr, "docker shim:", err)
		os.Exit(127)
	}
}
