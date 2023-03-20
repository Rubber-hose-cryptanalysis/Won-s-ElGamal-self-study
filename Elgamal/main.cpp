#include "ElGamal.hpp"

int main()
{
	std::vector<mpz_class> v;
	v.resize(4);

	keyGen(v);

	mpz_class m;
	mpz_set_str(m.get_mpz_t(), "1634976dbdeaba97a234f", 16);
	gmp_printf("The plaintext is %ZX\n", m);

	std::vector<mpz_class> c;
	c = Encrypt(m, v[0], v[1], v[2]);

	gmp_printf("The ciphertexts are (%ZX, %ZX)\n", c[0], c[1]);

	mpz_set(m.get_mpz_t(), Decrypt(c[0], c[1], v[3], v[2]).get_mpz_t());
	gmp_printf("The decrypted message is %ZX\n", m);

	return 0;
}