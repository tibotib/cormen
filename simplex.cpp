#include <bits/stdc++.h>
using namespace std;

template <typename T>
void display(vector<vector<T>> &A) {
        for(int i = 0; i < A.size(); ++i) {
                for(int j = 0; j < A[0].size(); ++j) cout << A[i][j] << ' ';
                cout << endl;
        }
        cout << endl;
}

template <typename T>
void display(vector<T> &vec) {
        for(int i = 0; i < vec.size(); ++i) cout << vec[i] << ' ';
        cout << endl<<endl;
}

template <typename T>
void pivot(vector<vector<T>> &A, vector<T> &b, vector<T> &c, T &nu,int e, int l) {
        A[l][e] = 1/A[l][e];
        b[l]    = -b[l] * A[l][e];
        for(int i = 0; i < A[0].size(); ++i) if(i != e) A[l][i] *= -A[l][e];

        for(int i = 0; i < A.size(); ++i) {
                int line = i;
                if(line == l) continue;

                for(int j = 0; j < A[0].size(); ++j) if(j != e) A[i][j] += A[i][e]*A[l][j];
                b[i]    += A[i][e] * b[l];
                A[i][e]  *= A[l][e];

        }

        for(int i = 0; i < A[0].size(); ++i) if(i != e) c[i] += A[l][i] * c[e];
        nu   += c[e]*b[l];
        c[e] *= A[l][e];
}

template <typename T>
int continue_simplex(const vector<T> &c) {
        for(int i = 0; i < c.size(); ++i) {
                if(c[i] > 0) return i;
        }
        return -1;
}

template <typename T>
vector<T> simplex(vector<vector<T>> A, vector<T> b, vector<T> c) {
        display(A);
        for(int i = 0; i < A.size(); ++i) for(int j = 0; j < A[0].size(); ++j) A[i][j] = -A[i][j];

        T nu = 0;
        int index = continue_simplex(c);
        while(index >= 0) {
                int eq = 0;
                T eq_value = (A[0][index] < 0) ?  -(b[0]/A[0][index]) : numeric_limits<T>::max();
                for(int i = 1; i < A.size(); ++i) {
                        if(A[i][index] >= 0) continue;
                        auto value = -(b[i] / A[i][index]);
                        if(value < eq_value) {
                                eq_value = value;
                                eq       = i;
                        }
                }
                if(eq_value == numeric_limits<T>::max()) {
                        cout <<"NON BORNE"<<endl;
                        break;
                }
                pivot(A,b,c,nu,index,eq);
                index = continue_simplex(c);
        }
        display(b);
        cout << "nu = " << nu << endl;

        return b;
}

int main() {/*
        vector<vector<double>> A{{-1,-1},{-1,-2},{3,5}};
        vector<double> c{1,1};
        vector<double> b{-5,-3,20};*/
        vector<vector<double>> A{{15,20,25, 59, 57, 78},{35,60,60,35,60,60},{20,30,25,35,60,60},{0,250,0, 0,250,0}};
        vector<double> c{300,250,450, 300, 250, 450};
        vector<double> b{1200,3000,1500,500};
        simplex<double>(A,b,c);
        return 0;
}
