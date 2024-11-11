#include <iostream>
#include <vector>

using namespace std;

class Matrix {
public: unsigned int l{}, w{};
    vector<vector<double>> v;

    Matrix() = default;

    Matrix(const unsigned int m, const unsigned int n) {
        this->l = m;
        this->w = n;
        this->v = vector<vector<double>>(m, vector<double>(n));
    }

    explicit Matrix(const vector<vector<double>>& t) {
        this->l = t.size();
        this->w = t[0].size();
        this->v = t;
    }

    Matrix operator+(const Matrix& k) const {
        if (k.l != this->l || k.w != this->w) {
            throw invalid_argument("Invalid dimension.");
        }

        Matrix m(this->l, k.w);
        for (int i = 0; i < this->l; i++) {
            for (int j = 0; j < k.w; j++) {
                m.v[i][j] = this->v[i][j] + k.v[i][j];
            }
        }
        return m;
    }

    bool operator==(const Matrix& k) const {
        if (this->l != k.l || this->w != k.w) return false;
        for (int i = 0; i < this->l; i++) {
            for (int j = 0; j < this->w; j++)
                if (this->v[i][j] != k.v[i][j]) return false;
        }
        return true;
    }

    bool operator!=(const Matrix& k) const {
        return !this->operator==(k);
    }

    Matrix operator*(const Matrix& k) const {
        if (this->w != k.l) {
            throw invalid_argument("Invalid dimension.");
        }

        double s = 0;

        Matrix m(this->l, k.w);
        for (int i = 0; i < this->l; i++) {
            for (int j = 0; j < k.w; j++) {
                for (int u = 0; u < this->w; u++) {
                    m.v[i][j] += this->v[i][u] * k.v[u][j];
                }
            }
        }
        return m;
    }

    vector<double> getDets(const Matrix& x) const {
        vector<double> d;
        for (int i = 0; i < this->w; i++) {
            d.push_back(this->replaceColWCol(x, i).det());
        }
        return d;
    }

    Matrix replaceColWCol(const Matrix& b, const int x) const {
        if (b.w != 1) {
            throw invalid_argument("Invalid dimension");
        }
        Matrix m(this->v);
        for (int i = 0; i < m.l; i++) {
            m.v[i][x] = b.v[i][x];
        }
        return m;
    }

    Matrix operator*(const double k) const {
        Matrix m(this->v);
        for (int i = 0; i < m.l; i++) {
            for (int j = 0; j < m.w; j++) {
                m.v[i][j] *= k;
            }
        }
        return m;
    }

    Matrix operator/(const double k) const {
        Matrix m(this->v);
        for (int i = 0; i < m.l; i++) {
            for (int j = 0; j < m.w; j++) {
                m.v[i][j] /= k;
            }
        }
        return m;
    }

    Matrix transpose() const {
        Matrix m(vector<vector<double>>(this->w, vector<double>(this->l)));
        for (int i = 0; i < m.w; i++) {
            for (int j = 0; j < m.l; j++) {
                m.v[i][j] = this->v[j][i];
            }
        }
        return m;
    }

    Matrix removeColAndRow(const int x, const int y) const {
        Matrix m(this->v);
        for (auto & i : m.v) {
            i.erase(i.begin() + y);
        }
        m.v.erase(m.v.begin() + x);
        m.l--;
        m.w--;
        return m;
    }

    static int getSign(const int i, const int j) {
        if ((i + j) % 2 == 0) return 1;
        return -1;
    }

    Matrix comp() const {
        Matrix m(this->l, this->w);
        for (int i = 0; i < m.l; i++) {
            for (int j = 0; j < m.w; j++) {
                m.v[i][j] = getSign(i, j) * this->removeColAndRow(i, j).det();
            }
        }
        return m;
    }

    double det() const {
        if (this->l != this->w) {
            throw invalid_argument("Given matrix is not square");
        }

        if (this->l == 1) return this->v[0][0];
        if (this->l == 2) return (this->v[0][0] * this->v[1][1]) - (this->v[0][1] * this->v[1][0]);

        double s = 0, p1, p2;

        vector<vector<double>> a(this->v);
        int n = a.size();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n - 1; j++) a[i].push_back(a[i][j]);
        }


        for (int k = 0; k < n; k++) {
            p1 = p2 = 1;
            for (int i = 0, j = k; i < (n + k) && j < (n + k); i++, j++) {
                p1 *= a[i][j];
            }

            for (int i = n - 1, j = k; i >= 0 && j < n + k; i--, j++) {
                p2 *= a[i][j];
            }
            s += p1 - p2;
        }
        return s;
    }

    Matrix adj() const {
        return this->comp().transpose();
    }

    Matrix inverse() const {
        if (this->l != this->w) {
            throw invalid_argument("Given matrix is not square");
        }
        if (this->det() == 0)
            throw invalid_argument("Inverse matrix non-existent");
        return this->adj() / this->det();
    }
};

int n;
Matrix a, b, s;
double d;

int main() {
    cin >> n;
    a = Matrix(n, n);
    b = Matrix(n, 1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) {
            if (j == n) cin >> b.v[i][0];
            else cin >> a.v[i][j];
        }
    }

    d = a.det();

    if (d == 0) {
        vector<double> dets = a.getDets(b);
        for (int i = 0; i < (int) dets.size(); i++) {
            if (dets[i] != 0) {
                cout << "No solution";
                return 0;
            }
        }
        cout << "The system of linear equations has infinitely many solutions:\n";

        return 0;
    }

    s = a.inverse() * (b);

    cout << "The system of linear equations has the only solution: \n";

    for (int i = 0; i < s.l; i++) {
        for (int j = 0; j < s.w; j++) cout << s.v[i][j] << ' ';
        cout << '\n';
    }

    return 0;
}