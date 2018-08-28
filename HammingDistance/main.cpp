// Calculate the count of "1" of an unsigned integer in binary expression.

#include <bitset>
#include <iostream>
using namespace std;

// method 01: count if the right bit of n + 1 is 1 until n is 0, then n = n/2 (bit right move)
int method01(int n) {
    int count(0);
    while (n != 0) {
        count += n & 1;
        n >>= 1;  // bit right move
    }
    return count;
}

// method 02 : n & (n-1) always equal zero, if large bit of n is 0, this method will be best
// Extension: How to check n is the power of 2 ? n > 0 && ((n & (n-1)) == 0
int method02(int n) {
    int count(0);
    while (n != 0) {
        n &= n - 1;
        ++count;
    }
    return count;
}

// method 03 : using the divide and conquer algorithm, please count the near 2 bits, then count the near 4 bits, and
// then 8 bits, 16 bits, 32 bits, refer to wiki
int method03(int n) {
    n = (n & 0x55555555) + ((n >> 1) & 0x55555555);
    n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
    n = (n & 0x0F0F0F0F) + ((n >> 4) & 0x0F0F0F0F);
    n = (n & 0x00FF00FF) + ((n >> 8) & 0x00FF00FF);
    n = (n & 0x0000FFFF) + ((n >> 16) & 0x0000FFFF);
    return n;
}

int main(int argc, char* argv[]) {
    int n(0);
    cout << "Please Input the number:";
    cin >> n;
    bitset<32> myBit(n);
    cout << "Number " << n << " = " << myBit << endl;
    cout << "Count of 1 in Binary Format: " << myBit.count() << endl;

    cout << "Method 01: " << method01(n) << endl;
    cout << "Method 02: " << method02(n) << endl;
    cout << "Method 03: " << method03(n) << endl;

    return 0;
}