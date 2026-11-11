#!/usr/bin/env bash
# BandQPI-3D v1.0 — build and check the installation
#   bash install.sh            build siesta3d.x here (gfortran + OpenMP + LAPACK)
#   bash install.sh --test     build, then run the synthetic-model test suite (exact answers)
set -e
cd "$(dirname "$0")"
FC=${FC:-gfortran}
command -v "$FC" >/dev/null || { echo "gfortran not found: install it (e.g. apt install gfortran, brew install gcc, conda install gfortran)"; exit 1; }
command -v python3 >/dev/null || { echo "python3 not found"; exit 1; }

# LAPACK: OpenBLAS if present, else reference LAPACK + BLAS, else MKL via LAPACK_LIBS
if [ -n "$LAPACK_LIBS" ]; then LIBS="$LAPACK_LIBS"
elif echo 'program t; end' > .t.f90 && $FC .t.f90 -lopenblas -o .t 2>/dev/null; then LIBS="-lopenblas"
elif $FC .t.f90 -llapack -lblas -o .t 2>/dev/null; then LIBS="-llapack -lblas"
else echo "no LAPACK found: install libopenblas-dev (or liblapack-dev), or set LAPACK_LIBS"; rm -f .t .t.f90; exit 1
fi
rm -f .t .t.f90
ARCH=${ARCH:--march=native}
echo "building siesta3d.x with $FC $LIBS"
$FC -O3 $ARCH -fopenmp -ffree-line-length-none siesta3d.f90 $LIBS -o siesta3d.x
rm -f *.mod

python3 - <<'EOF'
import importlib.util, sys
need = ['numpy', 'scipy', 'matplotlib', 'PIL']
miss = [m for m in need if importlib.util.find_spec(m) is None]
if miss:
    print('missing Python packages:', ' '.join(miss), '\n  pip install numpy scipy matplotlib pillow')
    sys.exit(1)
print('Python packages: numpy scipy matplotlib pillow  ok')
EOF
echo "installed: $(pwd)/siesta3d.x  (keep bandqpi3d.py in the same folder)"

if [ "$1" == "--test" ]; then
  # exact-answer test suite: ten synthetic SIESTA HSX models (about 10-15 minutes on 2 cores)
  set +e
  rm -rf tests/run && mkdir -p tests/run
  python3 tests/make_models.py tests/run
  for d in tests/run/*/; do
    (cd "$d" && ../../../siesta3d.x *.fdf --fast --noshow > run.log 2>&1
     echo "$(basename $d): exit $? | $(tr '\r' '\n' < run.log | grep -a 'validation   \|independent audit *: [0-9]' | tr -s ' ' | tr '\n' '|')")
  done
  python3 tests/check_models.py tests/run > /dev/null; s1=$?
  tail -n 1 tests/run/model_check.txt
  python3 tests/negative_controls.py tests/run/B_honeycomb_gap/honeycomb_gap_3dbands bandqpi3d.py; s2=$?
  if [ $s1 -eq 0 ] && [ $s2 -eq 0 ]; then echo "TEST SUITE PASSED (details: tests/run/model_check.txt)"; else echo "TEST SUITE FAILED"; exit 1; fi
fi
