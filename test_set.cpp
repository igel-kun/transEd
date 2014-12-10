#include <set>
#include <algorithm>
#include <vector>
#include <iostream>
using namespace std;

int main(){
	set<uint> s1, s2, result;
	set<uint>::iterator i1 = s1.begin();
	set<uint>::iterator i2 = s1.end();
	set<uint>::iterator i3 = s2.begin();
	set<uint>::iterator i4 = s2.end();
	set<uint>::iterator i5 = result.begin();
//	set<uint>::iterator it = set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), result.begin());
//	set<uint>::iterator it = set_union(i1, i2, i3, i4, i5);


	int first[] = {5,10,15,20,25};
	int second[] = {50,40,30,20,10};
	set<int> v;                           // 0  0  0  0  0  0  0  0  0  0
	set<int>::iterator it;
	
	sort (first,first+5);     //  5 10 15 20 25
	sort (second,second+5);   // 10 20 30 40 50
	
	it=set_union (first, first+5, second, second+5, v.begin());
	// 5 10 15 20 25 30 40 50  0  0
	cout << "union has " << v.size() << " elements.\n";
	return 0;
}
