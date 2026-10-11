# SelfAssembly
   Mente Carlo simulation of Self assembly of chiral molecular in the metal surface.

## Source variants

The three historical source files contain C++ despite having no extension:

| File | Variant |
| --- | --- |
| `chrial molecular +linker` | Chiral molecules with V-shaped linkers |
| `chrial molecular+C3 linker` | C3-symmetric molecules and annealing |
| `metal atom+linker` | Metal atoms and 120-degree linkers |

Compile a selected file explicitly as C++. For a GCC-compatible compiler:

```sh
g++ -x c++ -O2 "chrial molecular +linker" -o chiral-linker
```

Run each variant in its own output directory. The programs write relative
text outputs (`ndata*.txt`, `nnwdata*.txt`, or `data*.txt` depending on
the variant) and can replace existing files with the same names. L, particle
counts, schedules, and energy rules are constants in the source. Randomness
is seeded from the current time.

The chiral variants now initialize only valid array columns and wrap the
initial neighbor's column into [0,L). The scientific move rules are unchanged.
These repairs are not a complete bounds audit or a validation of simulation
results; review all lattice accesses before production simulations.

