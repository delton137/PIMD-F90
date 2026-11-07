# Regression test for the contracted-PIMD / SIESTA double-counting fix

`fsiesta_mock.f90` replaces `fsiesta_pipes.f90` with a mock of the SIESTA pipe
interface. It returns, for a set of non-interacting water molecules, the
Partridge-Schwenke (PS) monomer energy and forces in eV and eV/Ang, which is
what the monomer-corrected SIESTA (`E_DFT + E_PS - E_fit`, with `E_DFT ~ E_fit`)
returns for isolated molecules. Like the real SIESTA it receives the raw
coordinates, and it counts how many molecules arrive split across the box edge
(`mock_siesta_split_count.dat`: calls, molecules seen, molecules split).

Build: in the Makefile replace `fsiesta_pipes` by `fsiesta_mock` in `FILES`
(the monomer PES routines `pot_ttm/nasa.f90`, `nasa_mod.f90` and the other TTM
sources must be present) and run `make serial`. Then from a directory holding
`input_coords/`:

    ./PIMD_serial.x s_b1.inp     # 1 H2O, 1 bead,  contracted, 100 ps
    ./PIMD_serial.x s_b32.inp    # 1 H2O, 32 beads, contracted, 100 ps
    ./PIMD_serial.x pbc5.inp     # 5 H2O, 1 bead, 15 A box: box-crossing test

Expected H-DOS (`out_*_DOS.dat`) for one molecule at 300 K:

| code                    | bend      | stretch   | split molecules at SIESTA |
|-------------------------|-----------|-----------|---------------------------|
| before fix (Nov 2026)   | ~2230-2290 (= sqrt2 x 1611) | ~5550 (= sqrt2 x 3922) | 15-25% |
| after fix               | ~1550-1680 | ~3850-3950 | 0 |

With 32 beads the unfixed code shows the sqrt(2) "ghost" peaks near 2290 and
5400 cm^-1; the fixed code shows none.
