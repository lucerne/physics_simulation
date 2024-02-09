#include <unordered_set>
#include <tuple>

using namespace std;

int main () {
    unordered_set< tuple<int, int> > s ;
    tuple<int, int> t ;
    s.insert(t) ;
}