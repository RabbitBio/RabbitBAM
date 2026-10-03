#!/usr/bin/env python3
"""Build the internal MPI regression with an existing Unix Makefiles app build."""
import pathlib
import shlex
import subprocess
import sys

root = pathlib.Path(__file__).resolve().parents[1]
build = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else root / "build_sunway_sduhpc"
algorithm = sys.argv[2] if len(sys.argv) > 2 else "collate"
implementations = {
    "collate": ("collate_exchange_regression", "algorithms/collate/swbam_collate_mpi.cpp.o"),
    "sort": ("sort_exchange_regression", "algorithms/sort/swbam_sort_optimized_mpi.cpp.o"),
    "markdup": ("markdup_boundary_regression", "algorithms/markdup/swbam_markdup_mpi.cpp.o"),
}
if algorithm not in implementations:
    raise SystemExit("algorithm must be collate, sort, or markdup")
app = build / "swbam-app"
target = app / "CMakeFiles/RabbitBAM-MPI.dir"
name, implementation = implementations[algorithm]
output = build / name
settings = {}
for line in (target / "flags.make").read_text().splitlines():
    key, sep, value = line.partition(" = ")
    if sep:
        settings[key] = shlex.split(value)
link = shlex.split((target / "link.txt").read_text())
obj = str(output) + ".o"
subprocess.check_call([link[0]] + settings["CXX_DEFINES"] +
                      settings["CXX_INCLUDES"] + settings["CXX_FLAGS"] +
                      ["-c", str(root / f"tools/{name}.cpp"), "-o", obj],
                      cwd=app)
link[link.index("CMakeFiles/RabbitBAM-MPI.dir/main.cpp.o")] = obj
link.remove(f"CMakeFiles/RabbitBAM-MPI.dir/{implementation}")
link[link.index("-o") + 1] = str(output)
subprocess.check_call(link, cwd=app)
print(output)
