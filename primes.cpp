#include "CS207/Util.hpp"
#include <vector>
#include <algorithm> 
#include <math.h>

using namespace std;

/** Return true iff @a n is prime.
 * @pre @a n >= 0 */
bool is_prime(int n)
{
	int i;
  	assert(n >= 0);
 	
	static vector<int> v_primes;
	vector <int>::iterator Iter;

	// 	avoid segmentation fault using function "find" if container is empty
	if ( v_primes.empty() ) v_primes.push_back(2);

	Iter = find (v_primes.begin(), v_primes.end(), n);	
	
	
	// a good test is 5
	/* if n is the last element of the vector */
	if ( n == v_primes.back() ) return true;	
	/* if n is less than the last element of the vector, and iterator reaches the end of 
		of the vector */
		
	if ( n < v_primes.back() && Iter == v_primes.end() ) return false;
	/* if n is less than the last element of the vector, and iterator found something in 
		the vector */	
	if ( n < v_primes.back() && Iter != v_primes.end() ) return true;	
	
	/* if n is great than the last element of the vector, we need to check the integers
		past n */	
	if ( n > v_primes.back() ) {
	
	  	i = floor(sqrt(n));
  		while(n % i != 0 ){
  			i--;
  		}
  		
  		if (i == 1){
  			v_primes.push_back(n);
  			return true;
  		}
		else return false;	
		}
	
	}


int main()
{
  while (!cin.eof()) {
    // How many primes to test? And should we print them?
    cerr << "Input Number: ";
    int n = 0;
    CS207::getline_parsed(cin, n);
    if (n <= 0)
      break;

    cerr << "Print Primes (y/n): ";
    char confirm = 'n';
    CS207::getline_parsed(cin, confirm);
    bool print_primes = (confirm == 'y' || confirm == 'Y');

    CS207::Clock timer;

    // Loop and count primes from 2 up to n
    int num_primes = 0;
    for (int i = 2; i <= n; ++i) {
      if (is_prime(i)) {
        ++num_primes;
        if (print_primes)
          cout << i << endl;
      }
    }

    double elapsed_time = timer.elapsed();
    
    cout << "There are " << num_primes
         << " primes less than or equal to " << n << ".\n"
         << "Found in " << (1000 * elapsed_time) << " milliseconds.\n\n";
  }

  return 0;
}
