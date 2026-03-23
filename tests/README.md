# CLM-ml Functional Test Suite

Functional tests for individual Fortran subroutines and functions in CLM-ml v2.
Each test compiles a thin Fortran driver program that reads inputs from a
namelist (stdin) and writes `KEY=VALUE` lines to stdout.  A Python/pytest
layer runs the drivers with varied inputs and asserts on the outputs.

---

## Quick start

```bash
# 1. Build the main model (creates .o and .mod files in offline_executable/)
cd offline_executable
make

# 2. Build all test executables (links against the main-build objects)
cd ../tests
make all

# 3. Run the full test suite
cd python
python -m pytest -v

# 4. Run a single module's tests
python -m pytest test_water_vapor.py -v

# 5. Run tests matching a keyword
python -m pytest -k "satvap or latvap" -v
```

---

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| `nvfortran` on `PATH` | Same compiler as the main model |
| `NETCDF_PATH` set | e.g. `export NETCDF_PATH=/path/to/netcdf` |
| Main model built | Run `make` in `offline_executable/` first |
| Python ≥ 3.8 | Standard library only (no extra packages needed beyond pytest) |
| `pytest` installed | `pip install pytest` |

---

## Repository structure

```
tests/
├── Makefile                   # Builds test executables into tests/bin/
├── README.md                  # This file
├── .gitignore                 # Ignores bin/ and Python cache
├── fortran/                   # Fortran test-driver source files
│   ├── test_quadratic.F90
│   ├── test_tridiag.F90
│   ├── test_log_gamma.F90
│   ├── test_beta_function.F90
│   ├── test_beta_pdf.F90
│   ├── test_beta_cdf.F90
│   ├── test_SatVap.F90
│   ├── test_LatVap.F90
│   ├── test_ft.F90
│   ├── test_fth.F90
│   ├── test_fth25.F90
│   ├── test_phim_mo.F90       # phim_monin_obukhov + phic_monin_obukhov
│   ├── test_psim_mo.F90       # psim_monin_obukhov + psic_monin_obukhov
│   ├── test_shr_orb_params.F90
│   └── test_shr_orb_decl.F90
├── bin/                       # Compiled executables (created by make; git-ignored)
└── python/
    ├── conftest.py            # pytest session fixture: auto-runs make
    ├── utils.py               # run_fortran() / run_fortran_array() helpers
    ├── test_math_tools.py     # MLMathToolsMod tests
    ├── test_water_vapor.py    # MLWaterVaporMod tests
    ├── test_leaf_photo_temp.py# MLLeafPhotosynthesisMod: ft, fth, fth25
    ├── test_turbulence.py     # MLCanopyTurbulenceMod: Monin-Obukhov functions
    └── test_orb.py            # shr_orb_mod tests
```

---

## Covered subroutines / functions

| File | Module | Function/Subroutine | Test file |
|------|--------|---------------------|-----------|
| MLMathToolsMod.F90 | MLMathToolsMod | `quadratic` | test_math_tools.py |
| MLMathToolsMod.F90 | MLMathToolsMod | `tridiag` | test_math_tools.py |
| MLMathToolsMod.F90 | MLMathToolsMod | `log_gamma_function` | test_math_tools.py |
| MLMathToolsMod.F90 | MLMathToolsMod | `beta_function` | test_math_tools.py |
| MLMathToolsMod.F90 | MLMathToolsMod | `beta_distribution_pdf` | test_math_tools.py |
| MLMathToolsMod.F90 | MLMathToolsMod | `beta_distribution_cdf` | test_math_tools.py |
| MLWaterVaporMod.F90 | MLWaterVaporMod | `SatVap` | test_water_vapor.py |
| MLWaterVaporMod.F90 | MLWaterVaporMod | `LatVap` | test_water_vapor.py |
| MLLeafPhotosynthesisMod.F90 | MLLeafPhotosynthesisMod | `ft` | test_leaf_photo_temp.py |
| MLLeafPhotosynthesisMod.F90 | MLLeafPhotosynthesisMod | `fth` | test_leaf_photo_temp.py |
| MLLeafPhotosynthesisMod.F90 | MLLeafPhotosynthesisMod | `fth25` | test_leaf_photo_temp.py |
| MLCanopyTurbulenceMod.F90 | MLCanopyTurbulenceMod | `phim_monin_obukhov` | test_turbulence.py |
| MLCanopyTurbulenceMod.F90 | MLCanopyTurbulenceMod | `phic_monin_obukhov` | test_turbulence.py |
| MLCanopyTurbulenceMod.F90 | MLCanopyTurbulenceMod | `psim_monin_obukhov` | test_turbulence.py |
| MLCanopyTurbulenceMod.F90 | MLCanopyTurbulenceMod | `psic_monin_obukhov` | test_turbulence.py |
| shr_orb_mod.F90 | shr_orb_mod | `shr_orb_params` | test_orb.py |
| shr_orb_mod.F90 | shr_orb_mod | `shr_orb_decl` | test_orb.py |

> **Note on access visibility**: `ft`, `fth`, `fth25` (MLLeafPhotosynthesisMod) and
> the four Monin-Obukhov functions (MLCanopyTurbulenceMod) were originally
> `private`.  They have been changed to `public` in the source so they can be
> called from test drivers.  This is the only modification made to existing
> source files.

---

## How test drivers work

Each Fortran driver follows this pattern:

```fortran
program test_<name>
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog          ! CLM log unit
  use <Module>, only : <subroutine>
  implicit none
  real(r8) :: input1, output1
  namelist /inputs/ input1

  iulog = 6                             ! Direct CLM warnings to stdout
  input1 = <default>
  read(*, nml=inputs)                   ! Read inputs from stdin

  call <subroutine>(input1, output1)

  write(*, '(A,ES24.16)') 'output1=', output1   ! KEY=VALUE output
end program
```

Running manually:

```bash
echo "&inputs  t = 300.0  /" | ./bin/test_SatVap.exe
# Output:
# es=  3.5368...E+03
# desdt=  1.899...E+02
```

---

## How Python tests work

`utils.py` provides `run_fortran(exe_name, inputs_dict) -> dict`:

```python
from utils import run_fortran

out = run_fortran('test_SatVap.exe', {'t': 273.15})
print(out['es'])     # 611.2 Pa
print(out['desdt'])  # 44.4 Pa/K
```

For subroutines with array outputs (e.g. tridiag):

```python
from utils import run_fortran_array

u = run_fortran_array('test_tridiag.exe',
                      {'n': 3, 'a': [0.0,-1.0,-1.0], 'b': [2.0,2.0,2.0],
                       'c': [-1.0,-1.0,0.0], 'r': [1.0,0.0,1.0]},
                      base_key='u')
# u = [1.0, 1.0, 1.0]
```

---

## Adding a new test

### 1. Write the Fortran driver (`tests/fortran/test_<name>.F90`)

```fortran
program test_MyFunc
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MyModule, only : MyFunc
  implicit none
  real(r8) :: x, y
  namelist /inputs/ x
  iulog = 6
  x = 0._r8
  read(*, nml=inputs)
  call MyFunc(x, y)
  write(*, '(A,ES24.16)') 'y=', y
end program test_MyFunc
```

### 2. Add the target to the test `Makefile`

In `tests/Makefile`, append to `TARGETS`:

```makefile
TARGETS = \
  ...existing targets... \
  $(BIN)/test_MyFunc.exe
```

The pattern rule `$(BIN)/%.exe: $(SRC)/%.F90 $(MAIN_OBJS)` handles compilation
automatically — no extra rules needed for simple cases.

### 3. Write the Python test (`tests/python/test_<module>.py`)

```python
from utils import run_fortran

def test_my_func_reference():
    out = run_fortran('test_MyFunc.exe', {'x': 1.0})
    assert abs(out['y'] - 2.718) < 0.001   # Reference value

def test_my_func_property():
    y1 = run_fortran('test_MyFunc.exe', {'x': 1.0})['y']
    y2 = run_fortran('test_MyFunc.exe', {'x': 2.0})['y']
    assert y2 > y1   # Monotone increasing
```

---

## Troubleshooting

| Error | Likely cause | Fix |
|-------|-------------|-----|
| `FileNotFoundError: .exe not found` | Test not built yet | `make all` in `tests/` |
| `make: No rule to make target` | New `.F90` not in TARGETS | Add it to `tests/Makefile` |
| `Build failed: cannot find *.mod` | Main model not built | `make` in `offline_executable/` |
| `Build failed: NETCDF_PATH` | Environment not set | `export NETCDF_PATH=...` |
| Fortran `stop` / non-zero exit | `endrun` called in Fortran | Check stdout for error message |
| `SEGFAULT` or `NaN` | Uninitialised input | Check namelist defaults in driver |
