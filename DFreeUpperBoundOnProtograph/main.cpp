/*
Copyright(c) 2013, Ilya Vorobyev und Vasiliy Usatyuk
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met :
*Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and / or other materials provided with the distribution.
* Neither the name of the <organization> nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include<iostream>
#include<vector>
#include<queue>
#include<stack>
#include<set>
#include <iomanip>
#include <cassert>
#include<string>

using namespace std;
#define pii pair<int, int>
#define ll long long

float timeComb = 0;
float sumTime = 0;
float timePerm = 0;
bool nextCombination(vector<int>& a, int n) {
    int k = (int)a.size();
    for (int i = k - 1; i >= 0; --i)
        if (a[i] < n - k + i + 1) {
            ++a[i];
            for (int j = i + 1; j < k; ++j)
                a[j] = a[j - 1] + 1;
            return true;
        }
    return false;
}

//馬耶夫斯基混蛋

ll getPermanent(const vector<vector<int> >& sparse, const vector<vector<int> >& newMtr, vector<int>& usedRows, int er, int step = 0) {
    int n = sparse.size() - 1;
    ll res = 0;
    int col = step;
    if (col >= er)
        ++col;
    if (step + 1 == n) {
        for (int i = 0; i < sparse[col].size(); ++i) {
            if (usedRows[sparse[col][i]])
                continue;
            return newMtr[sparse[col][i]][col];
        }
    }
    for (int i = 0; i < sparse[col].size(); ++i) {
        if (usedRows[sparse[col][i]])
            continue;
        usedRows[sparse[col][i]] = 1;
        ll cur = getPermanent(sparse, newMtr, usedRows, er, step + 1);
        res += cur * newMtr[sparse[col][i]][col];
        usedRows[sparse[col][i]] = 0;

    }
    return res;
}





ll solve(const vector<int>& mask, const vector<vector<int> >& mtr, vector<vector<int> >& newMtr) {
    for (int i = 0; i < mtr.size(); ++i)
        for (int j = 0; j < mask.size(); ++j) {
            newMtr[i][j] = mtr[i][mask[j]];
        }
    ll res = 0;
    int n = newMtr.size();
    vector<vector<int> > sparseByCol(n + 1);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n + 1; ++j) {
            if (newMtr[i][j])
                sparseByCol[j].push_back(i);
        }
    }
    vector<int> usedRows(n, 0);
    
    for (int i = 0; i < mask.size(); ++i)
        res += getPermanent(sparseByCol, newMtr, usedRows, i);
    return res;
}



vector<vector<int> > sortByColWei(const vector<vector<int> >& mtr) {
    vector<int> colDeg(mtr[0].size(), 0);
    for (int i = 0; i < mtr.size(); ++i) {
        for (int j = 0; j < mtr[i].size(); ++j) {
            colDeg[j] += mtr[i][j];
        }
    }
    vector<pii> toSort;
    for (int i = 0; i < colDeg.size(); ++i) {
        toSort.push_back(pii(colDeg[i], i));
    }
    sort(toSort.begin(), toSort.end());
    vector<vector<int> > res(mtr.size(), vector<int>(mtr[0].size()));
    for (int i = 0; i < mtr.size(); ++i) {
        for (int j = 0; j < mtr[i].size(); ++j) {
            res[i][j] = mtr[i][toSort[j].second];
        }
    }
    return res;

}

ll countBound(vector<vector<int> > mtr) {
    mtr = sortByColWei(mtr);
    vector<vector<int> > newMtr(mtr.size(), vector<int>(mtr.size() + 1));
    int J = mtr.size(), I = mtr[0].size();
    if (I <= J)
        return -1;
    vector<int> mask(J + 1, 0);
    for (int i = 0; i < J + 1; ++i)
        mask[i] = i;
    ll res = -1;
    clock_t t1 = clock();
    do {
        clock_t t2 = clock();
        timeComb += (float)t2 - (float)t1;
        ll cur = solve(mask, mtr, newMtr);
        if (cur > 0) {
            if ((res < 0) || (cur < res))
                res = cur;
        }
        t1 = clock();
    } while (nextCombination(mask, I - 1));
    return res;
}

void print(const vector<int>& a) {
    for (int i = 0; i < a.size(); ++i)
        cout << a[i] << " ";
    cout << endl;
}


int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(0);
    string INPUT_FILENAME = "";
    string OUTPUT_FILENAME = "";
    for (int i = 1; i + 1 < argc; ++i) {
        if (string(argv[i]) == "-inputFile") {
            INPUT_FILENAME = argv[i + 1];
            ++i;
            continue;
        }
        if (string(argv[i]) == "-outputFile") {
            OUTPUT_FILENAME = argv[i + 1];
            ++i;
            continue;
        }
    }
    if (INPUT_FILENAME == "") {
        cerr << "wrong input\n";
        cerr << "Usage: DFreeUpperBoundOnProtograph.exe -inputFile INPUT.TXT -outputFile OUTPUT.TXT\n";
        return 0;
    }
    freopen(INPUT_FILENAME.c_str(), "r", stdin);

    int n, k;
    cin >> n >> k;
    vector<vector<int> > protograph(k, vector<int>(n));
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> protograph[i][j];
        }
    }
    clock_t t1 = clock();
    ll res8 = countBound(protograph);
    clock_t t2 = clock();
    sumTime = (float)t2 - float(t1);
    cerr << "sum time = " << sumTime / CLOCKS_PER_SEC << endl;
    cerr << "comb time = " << timeComb / CLOCKS_PER_SEC << endl;

    if (OUTPUT_FILENAME != "") {
        freopen(OUTPUT_FILENAME.c_str(), "a", stdout);
    }
    cout << INPUT_FILENAME << "\t" << res8 << endl;
    freopen(INPUT_FILENAME.c_str(), "a", stdout);
    cout << "\nVontobelTh8 Upper Bound = " << res8 << endl;
    return 0;
}
