#include <stdio.h>
#include  <time.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
int main()
{

//随机产生0~100的数。
     srand((int)time(NULL));
     int x=0;

    for (int i=1; i<=89; i++)
    {


        x=rand();
        cout<<x<<endl;

    }
    return 0;

}
