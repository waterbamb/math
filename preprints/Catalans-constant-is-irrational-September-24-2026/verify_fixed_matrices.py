"""Reproduce the fixed-matrix certificate in Section 4.4 using integer arithmetic modulo 101.

This script checks Table 1 only; it does not verify the full irrationality theorem."""
import json

def verify_certificate():
    P = 101
    inv = lambda x: pow(x % P, -1, P)
    m = [0] * 65
    m[0] = 2
    for i in range(2, 65, 2):
        m[i] = i * m[i - 2] * inv(i + 1) % P
    km = [0] * 65
    km[1] = 2
    kp = [0] * 59
    for d in range(2, 65):
        km[d] = ((d - 1) * km[d - 2] + m[d - 2] + m[d - 1]) * inv(d) % P
    for d in range(2, 59):
        kp[d] = ((d - 2) * kp[d - 2] + 2 * inv(d - 1)) * inv(d - 1) % P
    Y = [[0] * 59 for _ in range(65)]
    for i in range(65):
        Y[i][0] = km[i]
    for j in range(59):
        Y[0][j] = kp[j]
    for i in range(1, 65):
        for j in range(1, 59):
            Y[i][j] = (Y[i - 1][j - 1] - m[i - 1] * inv(j)) % P
    h1 = [0] * 65
    h2 = [0] * 65
    for i in range(1, 65):
        h1[i] = (h1[i - 1] + inv(i)) % P
        h2[i] = (h2[i - 1] + inv(i) ** 2) % P
    J = [[3 * inv(2) * (-h2[i] if i == j else (h1[i] - h1[j]) * inv(i - j)) % P for j in range(59)] for i in range(65)]
    E = [[0] * 65 for _ in range(45)]
    I = [[0] * 65 for _ in range(45)]
    E[0][62] = 1
    E[1][61] = 1
    I[1][62] = 1
    for a in (E, I):
        for d in range(2, 45):
            for k in range(65):
                a[d][k] = (2 * (a[d - 1][k + 1] if k < 64 else 0) - a[d - 2][k]) % P
    B = []
    for r in range(49):
        u = abs(r - 4)
        sgn = (r > 4) - (r < 4)
        y = [(E[u][k] - 2 * (E[u][k - 1] if k >= 1 else 0) + (E[u][k - 2] if k >= 2 else 0)) % P for k in range(65)]
        z = [sgn * (I[u][k] - 2 * (I[u][k - 1] if k >= 1 else 0) + (I[u][k - 2] if k >= 2 else 0)) % P for k in range(65)]
        v = [sum((y[i] * Y[i][j] - z[i] * J[i][j] for i in range(65))) % P for j in range(7, 59)]
        for _ in range(4):
            v = [(v[i] - v[i + 1]) % P for i in range(len(v) - 1)]
        B.append(v)
    expected = {0: [60, 68, 62, 79, 47, 32, 69, 57, 30, 72, 35, 66, 43, 20, 85, 48, 88, 4, 77, 54, 60, 79, 26, 68, 83, 39, 40, 65, 1, 68, 78, 24, 15, 98, 32, 22, 94, 9, 99, 10, 15, 75, 4, 2, 25, 53, 90, 79], 1: [38, 51, 90, 70, 5, 21, 88, 55, 45, 20, 35, 41, 77, 10, 18, 25, 76, 14, 38, 72, 6, 66, 56, 35, 83, 58, 56, 11, 17, 20, 30, 24, 28, 11, 25, 70, 79, 99, 66, 38, 4, 41, 63, 91, 63, 17, 90, 98], -1: [82, 92, 26, 87, 21, 74, 87, 88, 88, 3, 14, 23, 38, 58, 36, 20, 26, 33, 94, 74, 78, 45, 93, 86, 73, 66, 45, 30, 61, 3, 88, 27, 20, 58, 69, 48, 78, 39, 48, 1, 66, 18, 86, 93, 52, 92, 39, 77]}
    results = {}
    for sigma in (0, 1, -1):
        A = [[(B[r][k] + sigma * B[r + 1][k]) % P for k in range(48)] for r in range(48)]
        pivots = []
        swaps = []
        for i in range(48):
            if not A[i][i]:
                j = next((j for j in range(i + 1, 48) if A[j][i]))
                A[i], A[j] = (A[j], A[i])
                swaps.append([i, j])
            pivot = A[i][i]
            pivots.append(pivot)
            for j in range(i + 1, 48):
                scale = A[j][i] * inv(pivot) % P
                A[j] = [(a - scale * b) % P for a, b in zip(A[j], A[i])]
        results[str(sigma)] = {'pivots': pivots, 'swaps': swaps, 'matches_table': pivots == expected[sigma]}
        if pivots != expected[sigma]:
            raise ValueError(f'Pivot certificate mismatch for sigma={sigma}: {pivots}')
        if swaps != ([[30, 31], [45, 46]] if sigma == 1 else []):
            raise ValueError(f'Row-swap certificate mismatch for sigma={sigma}: {swaps}')
    return results
if __name__ == '__main__':
    print(json.dumps(verify_certificate(), indent=2))
    print('PASS: all 144 pivots and stated row swaps reproduced modulo 101')
