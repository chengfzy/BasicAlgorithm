//====================== Q01: Hamming Weight ==========================
//
// ����һ���ֽڵ��޷������ͱ�������������Ʊ�ʾ��1�ĸ���
//
// Application:
// 1. ͼƬ���ƶȼ��
// 2. ����ƥ��: Hamming Distance
//=====================================================================

#include <iostream>
#include <bitset>
using namespace std;

//����1:�����,����һ����������,��β��Ϊ1ʱ��1,Ȼ��/2(����),ֱ������Ϊ0Ϊֹ
int Method01(int n)
{
	int count(0);	//������������
	while (n != 0)
	{
		count += n & 1;
		n >>= 1;	//����
	}
	return count;
}

//����2: n��n-1��������λ��Զ��0, ֻ����1��λ��,��n�д������λ��0ʱ,���ø÷������
//���꣺�ж�һ�����Ƿ�Ϊ2���ݴη�, n > 0 && ((n & (n-1)) == 0
int Method02(int n)
{
	int count(0);
	while (n != 0)
	{
		n &= n - 1;
		++count;
	}
	return count;
}

//����3,Hamming Weight:ʹ�÷��ε�˼��,�ȼ���ÿ�����ڵ�2λ���м���1,�ټ���ÿ���ڵ�4λ���м���1,
//Ȼ��8λ,16λ,32λ
int Method03_HammingWeight(int n)
{
	n = (n & 0x55555555) + ((n >> 1) & 0x55555555);
	n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
	n = (n & 0x0F0F0F0F) + ((n >> 4) & 0x0F0F0F0F);
	n = (n & 0x00FF00FF) + ((n >> 8) & 0x00FF00FF);
	n = (n & 0x0000FFFF) + ((n >> 16) & 0x0000FFFF);
	return n;
}


int main()
{
	int n(0);
	cout << "Please Input the number:";
	cin >> n;
	bitset<32> myBit(n);
	cout << "Number " << n << " = " << myBit << endl;
	cout << "Count of 1 in Binary Format: " << myBit.count() << endl;

	cout << "Method 01: " << Method01(n) << endl;
	cout << "Method 02: " << Method02(n) << endl;
	cout << "Method 03, Hamming Weight: " << Method03_HammingWeight(n) << endl;

	system("pause");
	return 0;
}