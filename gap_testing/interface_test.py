from subprocess import Popen, PIPE
gap_command = "/Users/chris/Downloads/gap-4.11.1/bin/gap.sh"

process = Popen([gap_command], stdout=PIPE, stdin=PIPE)
process.stdin.write(b"2+2;")
process.stdin.flush()

out = process.stdout.read()
print(repr(out))
process.stdin.write(b'2*2;')
process.stdin.flush()

