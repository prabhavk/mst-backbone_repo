#pragma once
#include <cstddef> 
#include <initializer_list>
#include <algorithm>
#include <ostream>
#include <iomanip>

struct Matrix4f {
    float data[4][4];

    // Default constructor - does not zero
    Matrix4f() = default;

    // Initialize from value
    explicit Matrix4f(float val) {
        for (size_t i = 0; i < 4; ++i)
            for (size_t j = 0; j < 4; ++j)
                data[i][j] = val;
    }

    // Element access with (row, col)
    float& operator()(size_t row, size_t col) {
        return data[row][col];
    }

    const float& operator()(size_t row, size_t col) const {
        return data[row][col];
    }

    // Set all entries to zero
    void setZero() {
        for (size_t i = 0; i < 4; ++i)
            for (size_t j = 0; j < 4; ++j)
                data[i][j] = 0.0f;
    }

    // Static method like Eigen::ArrayXXf::Zero(4,4)
    static Matrix4f Zero() {
        Matrix4f m;
        m.setZero();
        return m;
    }
    
};

struct Matrix4f {
    // storage
    std::array<std::array<float, 4>, 4> a{};  // value-initialized to 0s

    // row access so you can do M[i][j]
    float* operator[](int i) { return a[static_cast<size_t>(i)].data(); }
    const float* operator[](int i) const { return a[static_cast<size_t>(i)].data(); }

    // helpers
    static Matrix4f Zero() { return Matrix4f{}; }

    static Matrix4f Identity() {
        Matrix4f M;
        for (int i = 0; i < 4; ++i) M[i][i] = 1.0f;
        return M;
    }
};

// transpose
inline Matrix4f transposeMatrix4f(const Matrix4f& m) {
    Matrix4f t;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            t[i][j] = m[j][i];
    return t;
}

// stream output
inline std::ostream& operator<<(std::ostream& os, const Matrix4f& M) {
    auto old_flags = os.flags();
    auto old_prec  = os.precision();
    os.setf(std::ios::fixed, std::ios::floatfield);
    os.precision(6);

    for (int i = 0; i < 4; ++i) {
        os << "[ ";
        for (int j = 0; j < 4; ++j) {
            os << M[i][j];
            if (j < 3) os << ' ';
        }
        os << " ]";
        if (i < 3) os << '\n';
    }

    os.flags(old_flags);
    os.precision(old_prec);
    return os;
}

inline Matrix4f operator*(const Matrix4f& A, const Matrix4f& B) {
    Matrix4f C; // value-initialized to zeros thanks to std::array
    for (int i = 0; i < 4; ++i) {
        for (int k = 0; k < 4; ++k) {
            const float aik = A[i][k];
            for (int j = 0; j < 4; ++j) {
                C[i][j] += aik * B[k][j];
            }
        }
    }
    return C;
}

inline Matrix4f& operator*=(Matrix4f& A, const Matrix4f& B) {
    A = A * B;
    return A;
}

// (optional) scalar multiply if you ever need it
inline Matrix4f operator*(const Matrix4f& A, float s) {
    Matrix4f C;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            C[i][j] = A[i][j] * s;
    return C;
}
inline Matrix4f operator*(float s, const Matrix4f& A) { return A * s; }