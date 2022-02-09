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
void pivot(vector<vector<T>> &A, vector<T> &b, vector<T> &c, T &nu, vector<int> &N, vector<int> &B, int ee, int ll, map<int, int> &mp) {
        int l = mp[ll];
        int e = mp[ee];
        A[l][e] = 1/A[l][e];
        b[l]    = -b[l] * A[l][e];
        for(int i = 0; i < N.size(); ++i) {
                if(mp[N[i]] == e) continue;
                A[l][mp[N[i]]] *= -A[l][e];
        }

        for(int i = 0; i < B.size(); ++i) {
                int line = mp[B[i]];
                if(line == l) continue;

                auto v = A[line][e];
                for(int j = 0; j < N.size(); ++j) {
                        if(mp[N[j]] == e) continue;
                        A[line][mp[N[j]]] += v*A[l][mp[N[j]]];
                }
                A[line][e] = v*A[l][e];
                b[line] += v * b[l];
        }

        auto v = c[e];
        for(int i = 0; i < N.size(); ++i) {
                if(mp[N[i]] == e) continue;
                c[mp[N[i]]] += A[l][mp[N[i]]] * v;
        }

        c[e] = A[l][e] * v;
        nu += v*b[l];

        for(int i = 0; i < N.size(); ++i) {
                if(N[i] == ee) N[i] = ll;
        }
        for(int i = 0; i < B.size(); ++i) if(B[i] == ll) B[i] = ee;
        int t = mp[ee];

        mp[ee] = mp[ll];
        mp[ll] = t;
}

template <typename T>
int continue_simplex(const vector<T> &c, const vector<int> &N, map<int,int>&mp) {
        for(int i = 0; i < N.size(); ++i) {
                if(c[mp[N[i]]] > 0) return N[i];
        }
        return -1;
}

template <typename T>
vector<T> simplex(vector<vector<T>> A, vector<T> b, vector<T> c) {
        map<int, int> mp;
        vector<int> N(A[0].size());
        vector<int> B(A.size());
        for(int i = 0; i < A[0].size();++i) {
                mp.insert({i, i});
                N[i] =  i;
        }
        for(int i = 0; i < A.size(); ++i) {
                mp.insert({i + A[0].size(), i});
                B[i] = i+A[0].size();
        }
        for(int i = 0; i < A.size(); ++i) {
                for(int j =  0; j<  A[0].size(); ++j) A[i][j] = -A[i][j];
        }
        T nu = 0;
        int index = continue_simplex(c,N,mp);
        while(index >= 0) {
                int eq = 0;
                T eq_value = (A[mp[B[0]]][mp[index]] < 0) ?  -(b[mp[B[0]]]/A[mp[B[0]]][mp[index]]) : numeric_limits<T>::max();
                for(int i = 1; i < B.size(); ++i) {
                        if(A[mp[B[i]]][mp[index]] >= 0) continue;
                        auto value = -(b[mp[B[i]]] / A[mp[B[i]]][mp[index]]);
                        if(value < eq_value) {
                                eq_value = value;
                                eq = B[i];
                        }
                }
                if(eq_value == numeric_limits<T>::max()) {
                        cout <<"NON BORNE"<<endl;
                        break;
                }
                pivot(A,b,c,nu,N,B,index,eq,mp);
                index = continue_simplex(c,N,mp);
        }
        return b;
}

int main() {
        /*vector<vector<double>> A{{1,1,3},{2,2,5},{4,1,2}};
        vector<double> c{3,1,2};
        vector<double> b{30,24,36};*/
        vector<vector<double>> A{{1,1,-1},{-1,-1,1},{1,-2,2}};
        vector<double> c{2,-3,3};
        vector<double> b{7,-7,4};
        simplex<double>(A,b,c);
        return 0;
}
