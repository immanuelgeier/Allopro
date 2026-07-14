import numpy as np

n = 5

# zufällige symmetrische 3x3-Matrizen
B = np.random.randn(n, n, 3, 3)

for k in range(n):
    for l in range(n):
        B[k, l] = 0.5 * (B[k, l] + B[k, l].T)

x = np.random.randn(n)

y, sigma = find_best_y(B, x)

A = build_A(B, x, y)

print("kleinster Singulärwert =", sigma)
print()
print("A =")
print(A)
print()
print("Off-Diagonale:")
print(A[0,1], A[0,2], A[1,2])

def find_best_y(B, x):
    """
    Parameters
    ----------
    B : ndarray, shape (n,n,3,3)
        B[k,l] ist eine 3x3-Matrix.
    x : ndarray, shape (n,)
        Vorgegebener Vektor.

    Returns
    -------
    y : ndarray, shape (n,)
        Normierter Vektor, der die Off-Diagonalelemente
        von A möglichst klein macht.
    sigma_min : float
        Kleinster Singulärwert.
    """

    n = len(x)

    # relevante Off-Diagonal-Indizes
    pairs = [(0, 1), (0, 2), (1, 2)]

    C = np.zeros((3, n))

    for row, (i, j) in enumerate(pairs):

        # M^(ij)
        M = B[:, :, i, j]

        # c^(ij) = M^T x
        C[row] = M.T @ x

    # SVD
    U, S, Vh = np.linalg.svd(C)

    # rechter Singulärvektor zum kleinsten Singulärwert
    y = Vh[-1]

    return y / np.linalg.norm(y), S[-1]


def build_A(B, x, y):
    """
    Berechnet
        A = sum_{k,l} x_k y_l B_{k,l}.
    """

    A = np.tensordot(
        np.outer(x, y),
        B,
        axes=([0, 1], [0, 1])
    )

    return A