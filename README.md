د# Linear Algebra Project in C++

## File Structure

```text
Linear-Algebra-Project/
│
├── README.md
├── solve_system.cpp
├── transpose.cpp
├── determinant.cpp
├── inverse.cpp
├── adjoint.cpp
├── cofactors.cpp
└── minors.cpp
```

---

# solve_system.cpp

```cpp
#include <iostream>
#include <cmath>
using namespace std;

int main() {
    const int N = 10;

    double A[N][N + 1] = {
        { 2, -1,  3, -4,  5,  1, -2,  3,  4, -1,  17},
        { 1,  3, -2,  4, -5, -2,  3, -4, -5,  2, -23},
        { 3, -2,  4,  5, -1, -3,  4, -5,  2, -4,  14},
        { 4,  5, -1, -3,  4, -5,  2, -4,  5, -2, -15},
        {-1,  4,  5,  4, -3,  2, -4,  5, -2,  3,  11},
        { 2, -5, -4,  3,  2, -5,  4, -3,  2,  5,  -7},
        {-3,  4, -5, -1,  4,  2, -5,  4,  2,  3,  12},
        { 5, -2,  4,  5, -3, -4,  5, -2,  4, -1,  -8},
        { 4,  3, -2, -4,  5, -1,  3, -2, -1,  4, -19},
        {-1,  2,  3, -4, -5,  2, -3,  4,  5, -2,  10}
    };

    for (int col = 0; col < N; col++) {
        if (fabs(A[col][col]) < 1e-12) {
            cout << "Pivot is zero." << endl;
            return 0;
        }

        for (int row = col + 1; row < N; row++) {
            double factor = A[row][col] / A[col][col];

            for (int k = col; k <= N; k++) {
                A[row][k] -= factor * A[col][k];
            }
        }
    }

    double x[N];

    for (int row = N - 1; row >= 0; row--) {
        double sum = A[row][N];

        for (int col = row + 1; col < N; col++) {
            sum -= A[row][col] * x[col];
        }

        x[row] = sum / A[row][row];
    }

    char name[N] = {'a','b','c','d','e','f','g','h','i','j'};

    for (int i = 0; i < N; i++) {
        cout << name[i] << " = " << x[i] << endl;
    }

    return 0;
}
```

---

# transpose.cpp

```cpp
#include <iostream>
using namespace std;

int main() {
    const int N = 7;

    int A[N][N] = {
        {6, 4, -8, 9, 7, 6, 5},
        {7, 8, 6, -5, 4, 7, 8},
        {4, 0, 7, 8, 9, 7, 6},
        {5, 6, 4, 5, 7, 8, 7},
        {5, 4, -2, -5, 7, 3, -2},
        {4, 5, 7, 8, 9, 7, 6},
        {3, 2, 4, 5, 7, 8, 6}
    };

    int AT[N][N];

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            AT[j][i] = A[i][j];
        }
    }

    cout << "Transpose Matrix:" << endl;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << AT[i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}
```

---

# determinant.cpp

```cpp
#include <iostream>
#include <cmath>
using namespace std;

int main() {
    const int N = 7;

    double A[N][N] = {
        {6, 4, -8, 9, 7, 6, 5},
        {7, 8, 6, -5, 4, 7, 8},
        {4, 0, 7, 8, 9, 7, 6},
        {5, 6, 4, 5, 7, 8, 7},
        {5, 4, -2, -5, 7, 3, -2},
        {4, 5, 7, 8, 9, 7, 6},
        {3, 2, 4, 5, 7, 8, 6}
    };

    double det = 1.0;

    for (int col = 0; col < N; col++) {

        if (fabs(A[col][col]) < 1e-12) {
            cout << "Pivot is zero." << endl;
            return 0;
        }

        for (int row = col + 1; row < N; row++) {
            double factor = A[row][col] / A[col][col];

            for (int k = col; k < N; k++) {
                A[row][k] -= factor * A[col][k];
            }
        }
    }

    for (int i = 0; i < N; i++) {
        det *= A[i][i];
    }

    cout << "det(A) = " << llround(det) << endl;

    return 0;
}
```

---

# inverse.cpp

```cpp
#include <iostream>
#include <cmath>
using namespace std;

int main() {
    const int N = 7;

    double A[N][2 * N] = {
        {6,4,-8,9,7,6,5,   1,0,0,0,0,0,0},
        {7,8,6,-5,4,7,8,  0,1,0,0,0,0,0},
        {4,0,7,8,9,7,6,   0,0,1,0,0,0,0},
        {5,6,4,5,7,8,7,   0,0,0,1,0,0,0},
        {5,4,-2,-5,7,3,-2,0,0,0,0,1,0,0},
        {4,5,7,8,9,7,6,   0,0,0,0,0,1,0},
        {3,2,4,5,7,8,6,   0,0,0,0,0,0,1}
    };

    for (int col = 0; col < N; col++) {

        if (fabs(A[col][col]) < 1e-12) {
            cout << "Matrix is singular." << endl;
            return 0;
        }

        double pivot = A[col][col];

        for (int j = 0; j < 2 * N; j++) {
            A[col][j] /= pivot;
        }

        for (int row = 0; row < N; row++) {

            if (row != col) {

                double factor = A[row][col];

                for (int j = 0; j < 2 * N; j++) {
                    A[row][j] -= factor * A[col][j];
                }
            }
        }
    }

    cout << "Inverse Matrix A^-1:" << endl;

    for (int i = 0; i < N; i++) {

        for (int j = N; j < 2 * N; j++) {
            cout << A[i][j] << " ";
        }

        cout << endl;
    }

    return 0;
}
```
