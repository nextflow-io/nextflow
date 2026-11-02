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

package main

import (
	"reflect"
	"testing"
)

func TestTranslateVolume(t *testing.T) {
	cases := map[string]string{
		`C:\work\ab:/c/work/ab`:           `/mnt/c/work/ab:/c/work/ab`,
		`D:\My Data:/d/My Data:ro,Z`:      `/mnt/d/My Data:/d/My Data:ro,Z`,
		`d:/data:/d/data`:                 `/mnt/d/data:/d/data`,
		`C:\:/c`:                          `/mnt/c/:/c`,
		`/home/user/work:/home/user/work`: `/home/user/work:/home/user/work`,
		`/tmp:/tmp:ro`:                    `/tmp:/tmp:ro`,
	}
	for spec, expected := range cases {
		if got := translateVolume(spec); got != expected {
			t.Errorf("translateVolume(%q) = %q, expected %q", spec, got, expected)
		}
	}
}

func TestTranslateArgs(t *testing.T) {
	args := []string{"run", "--rm", "-v", `C:\a:/c/a`, "--volume", `D:\b:/d/b:ro`, "--volume=E:\\c:/e/c", "-w", "/c/a", "ubuntu", "-v"}
	expected := []string{"run", "--rm", "-v", `/mnt/c/a:/c/a`, "--volume", `/mnt/d/b:/d/b:ro`, "--volume=/mnt/e/c:/e/c", "-w", "/c/a", "ubuntu", "-v"}
	if got := translateArgs(args); !reflect.DeepEqual(got, expected) {
		t.Errorf("translateArgs() = %q, expected %q", got, expected)
	}
}
