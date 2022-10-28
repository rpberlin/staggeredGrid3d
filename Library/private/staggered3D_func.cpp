#include <iostream>
#include <vector>
using namespace std;

void printSomething(void)
{
  cout<<"Hello World"<<endl;
  return;
}

void display(vector<vector<vector<char> > >& ch)
{
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                cout << "ch[" << i << "][" << j << "]["
                     << k << "] = " << ch[i][j][k] << endl;
            }
        }
    }
    return;
}
