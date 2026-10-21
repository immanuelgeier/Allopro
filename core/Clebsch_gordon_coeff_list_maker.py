#!/usr/bin/env python3

import sys
from sympy.physics.wigner import clebsch_gordan


def generate_cg_file(lmax, filename):
    """
    Generate all Clebsch-Gordan coefficients

        <l1 m1, l2 m2 | l3 m3>

    with l1, l2, l3 <= lmax.

    The coefficients are written as double precision floating point
    numbers for easy reading from Fortran.
    """

    count = 0

    with open(filename, "w") as f:

        # Header
        f.write("# Clebsch-Gordan coefficients\n")
        f.write(f"# lmax = {lmax}\n")
        f.write("# Format: l1 m1 l2 m2 l3 m3 CG\n")

        for l1 in range(lmax + 1):
            for l2 in range(lmax + 1):

                # Triangle condition
                l3_min = abs(l1 - l2)
                l3_max = min(l1 + l2, lmax)

                for l3 in range(l3_min, l3_max + 1):

                    for m1 in range(-l1, l1 + 1):
                        for m2 in range(-l2, l2 + 1):

                            m3 = m1 + m2

                            # m3 must be allowed
                            if abs(m3) > l3:
                                continue

                            cg = clebsch_gordan(
                                l1, l2, l3,
                                m1, m2, m3
                            )

                            # Convert SymPy exact expression to double
                            cg_float = float(cg.evalf(17))

                            f.write(
                                f"{l1:4d} {m1:4d} "
                                f"{l2:4d} {m2:4d} "
                                f"{l3:4d} {m3:4d} "
                                f"{cg_float:.17e}\n"
                            )

                            count += 1

    print(f"Wrote {count} coefficients to {filename}")


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage:")
        print("    python generate_cg.py LMAX output.txt")
        print()
        print("Example:")
        print("    python generate_cg.py 10 clebsch_gordan.txt")
        sys.exit(1)

    lmax = int(sys.argv[1])
    filename = sys.argv[2]

    generate_cg_file(lmax, filename)