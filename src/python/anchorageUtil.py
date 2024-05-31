"""
Description:
    Some utilities for anchorage.

BSD 3-Clause License

Copyright (c) 2024, Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
from sys import argv, stderr
import subprocess


def RunShellCommand(cmd: str, msg: str, ignoreError: bool = False) -> list[str]:
    if(msg != ""):
        print(msg)
        # print(msg, "Command: {}".format(cmd))

    execute = subprocess.run(cmd.split(), capture_output=True, text=True)
    # print("stdout:", execute.stdout)
    # print("stderr:", execute.stderr)

    if not ignoreError and execute.returncode != 0:
        print("Error encountered when executing " + msg , file=stderr)
        print(execute)
        print("stdout:", execute.stdout, file=stderr)
        print("stderr:", execute.stderr, file=stderr)
        exit(1)

    return execute.stdout


__REV_COMP_TABLE__ = str.maketrans("ACGTURYSWKMBDHVN", "TGCAAYRSWMKVHDBN")
def revcomp(seq: str) -> str:
    global __REV_COMP_TABLE__
    rc_seq: str = seq.translate(__REV_COMP_TABLE__)[::-1]
    return rc_seq


def __anchorageUtil_revcomp_test_():
    s = "ATCGATCGATCG"
    r = revcomp(s)
    assert r == "CGATCGATCGAT"
    print("s:",s, "\nr:", r)

    s = "AAAAAAAA"
    r = revcomp(s)
    assert r == "TTTTTTTT"
    print("s:",s, "\nr:", r)
    
    s = "UUUUUU"
    r = revcomp(s)
    assert r == "AAAAAA"
    print("s:",s, "\nr:", r)

    s = "TTTTTT"
    r = revcomp(s)
    assert r == "AAAAAA"
    print("s:",s, "\nr:", r)
    
    s = "BVHDCUASW"
    r = revcomp(s)
    assert r == "WSTAGHDBV"
    print("s:",s, "\nr:", r)

    print("Debug tests for anchorageUtil.revcomp completed")
    return 0


def __anchorageUtil_RunShellCommand_test_():
    RunShellCommand("echo ATCG GCTA 1234", "")
    print("Debug tests for anchorageUtil.RunShellCommand completed")
    return 0

def __anchorageUtil_test_():
    __anchorageUtil_revcomp_test_()
    __anchorageUtil_RunShellCommand_test_()
    print("Debug tests for anchorageUtil completed")
    return 0

if __name__ == "__main__":
    if len(argv) >= 2 and argv[1] in ['test', '--test', '-test']:
        print("Debug tests for anchorageUtil:")
        __anchorageUtil_test_()
    else:
        raise NotImplementedError("You should not run anchorageUtil from main. Import and call the object ")