#include "ElGamal.hpp"

std::vector<mpz_class> sieve_of_Eratosthenes(mpz_class n)
{
	mpz_class temp; mpz_t r; mpz_init(r); mpz_sqrtrem(temp.get_mpz_t(), r, n.get_mpz_t()); //temp=floor(sqrt(n))

	unsigned long long finish = 0;
	mpz_export(&finish, 0, -1, sizeof finish, 0, 0, temp.get_mpz_t()); //finish=temp (mpz to unsigned long long)
	
	std::vector<bool> toinit(finish + 1, true);
	std::vector<std::vector<bool>> tocheck(finish + 1, toinit);
	std::vector<mpz_class> primes;

	tocheck[0][0] = false; tocheck[0][1] = false; //0,1은 소수가 아님.

	mpz_class i;
	for (mpz_set_ui(i.get_mpz_t(), 2); mpz_cmp(i.get_mpz_t(), temp.get_mpz_t())<=0; i++) //2~floor(sqrt(n))까지
	{
		if (!tocheck[0][mpz_get_ui(i.get_mpz_t())])
			continue;

		mpz_class j;
		for (mpz_mul(j.get_mpz_t(), i.get_mpz_t(), i.get_mpz_t()); mpz_cmp(j.get_mpz_t(), n.get_mpz_t()) <= 0; mpz_add(j.get_mpz_t(), j.get_mpz_t(), i.get_mpz_t()))
		{
			mpz_t q; mpz_init(q);
			mpz_fdiv_qr_ui(q, r, j.get_mpz_t(), finish + 1);

			tocheck[mpz_get_ui(q)][mpz_get_ui(r)] = false;
		}

	}

	for (unsigned long long s = 0; s <= finish; s++)
	{
		for (unsigned long long t = 0; t <= finish; t++)
		{
			if (tocheck[s][t])
			{
				mpz_set_ui(temp.get_mpz_t(), finish + 1);
				mpz_mul_ui(temp.get_mpz_t(), temp.get_mpz_t(), s);
				mpz_add_ui(temp.get_mpz_t(), temp.get_mpz_t(), t);

				if (mpz_cmp(temp.get_mpz_t(), n.get_mpz_t()) > 0)
					return primes;

				primes.emplace_back(temp);
			}
		}
	}

	mpz_clear(r);
	return primes;
}

std::vector<mpz_class> primes_lowerthan(mpz_class n)
{
	std::vector<mpz_class> v;
	mpz_class p; mpz_set_ui(p.get_mpz_t(), 2); v.emplace_back(p);

	while (mpz_cmp(p.get_mpz_t(), n.get_mpz_t()) <= 0)
	{
		mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
		v.emplace_back(p);
	}

	return v;
}

mpz_class Pollard_Rho(mpz_class n)
{
	mpz_t x, y; mpz_init(x); mpz_init(y);
	mpz_class d;

	mpz_set_ui(x, 2); mpz_set_ui(y, 2); mpz_set_ui(d.get_mpz_t(), 1);

	while (mpz_cmp_ui(d.get_mpz_t(), 1) == 0)
	{
		mpz_mul(x, x, x); mpz_add_ui(x, x, 1); //f(x)=x^2+1, x<-f(x)
		mpz_mul(y, y, y); mpz_add_ui(y, y, 1); mpz_mul(y, y, y); mpz_add_ui(y, y, 1); //y<=f(f(y))

		mpz_sub(d.get_mpz_t(), x, y); mpz_abs(d.get_mpz_t(), d.get_mpz_t());
		mpz_gcd(d.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t()); //d<-gcd(|x-y|,n)
	}

	mpz_clear(x); mpz_clear(y);
	return d;
}

mpz_class primitive_root(mpz_class p)
{
	mpz_class a;
	mpz_set_ui(a.get_mpz_t(), 2);

	mpz_t q; mpz_init(q);
	mpz_sub_ui(q, p.get_mpz_t(), 1); mpz_divexact_ui(q, q, 2); //q=(p-1)/2.

	mpz_t e; mpz_init(e);
	int i = 1;
	while (i < 3)
	{
		if (i == 1)
		{
			mpz_powm_ui(e, a.get_mpz_t(), 2, p.get_mpz_t());

			if (mpz_cmp_ui(e, 1) == 0)
			{
				i = 0;
				mpz_add_ui(a.get_mpz_t(), a.get_mpz_t(), 1);
			}
		}

		else if (i == 2)
		{
			mpz_powm(e, a.get_mpz_t(), q, p.get_mpz_t());

			if (mpz_cmp_ui(e, 1) == 0)
			{
				i = 0;
				mpz_add_ui(a.get_mpz_t(), a.get_mpz_t(), 1);
			}
		}

		i++;
	}

	mpz_clear(q); mpz_clear(e);
	return a;
}

bool subsafe_check(mpz_class n)
{
	/*
	int bound = (int)pow(2, 16);
	int finish = (int)pow(2, 8);

	std::vector<bool> tocheck(bound + 1, true);
	tocheck[0] = false; tocheck[1] = false;

	for (int i = 2; i <= finish; i++)
	{
		if (!tocheck[i])
			continue;

		for (int j = i * i; j <= bound; j += i)
			tocheck[j] = false;
	}
	*/
	std::vector<unsigned long> primes;
	/*
	for (int i = 2; i <= bound; i++)
		if (tocheck[i])
			primes.emplace_back(i);
	//primes에 2^16까지의 모든 소수를 저장.
	*/

	FILE* fp = fopen("primes.txt", "r");
	if (!fp)
		return false;

	int temp;
	while (fscanf(fp, "%d", &temp) != EOF)
		primes.emplace_back(temp);

	fclose(fp);

	for (int i = 1; i < size(primes); i++)
	{
		if (mpz_congruent_ui_p(n.get_mpz_t(), 1, primes[i]) != 0)
			return false;
	}

	return true;
}


std::vector<mpz_class> keyGen(std::vector<mpz_class>& v)
{
	mpz_class g, h, p, d;
	mpz_t q; mpz_init(q);

	gmp_randstate_t state; gmp_randinit_mt(state); gmp_randseed_ui(state, (96889010407 * time(NULL) + 1013904223) % (unsigned long)pow(2, 31) - 1);
	mpz_t temp; mpz_init(temp);

	mpz_set_ui(p.get_mpz_t(), 1);
	for (int i = 0; i < 2047; i++)
	{
		mpz_mul_2exp(p.get_mpz_t(), p.get_mpz_t(), 1);

		mpz_urandomb(temp, state, 1);
		mpz_add(p.get_mpz_t(), p.get_mpz_t(), temp);
	}
	
	mpz_set_ui(q, 1);
	while (mpz_probab_prime_p(q, 50) == 0)
	{
		mpz_nextprime(p.get_mpz_t(), p.get_mpz_t()); //p는 2048bit prime
		while (!subsafe_check(p))
			mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
		mpz_sub_ui(q, p.get_mpz_t(), 1); mpz_divexact_ui(q, q, 2); //p=2*q+1. p가 safe prime이면 완료.
	}

	mpz_set_ui(d.get_mpz_t(), 1);
	for (int i = 0; i < 239; i++)
	{
		mpz_mul_2exp(d.get_mpz_t(), d.get_mpz_t(), 1);

		mpz_urandomb(temp, state, 1);
		mpz_add(d.get_mpz_t(), d.get_mpz_t(), temp);
	}
	//d는 240bit number
	
	mpz_set(g.get_mpz_t(), primitive_root(p).get_mpz_t()); //g는 p의 primitive root
	mpz_powm(h.get_mpz_t(), g.get_mpz_t(), d.get_mpz_t(), p.get_mpz_t()); //h = g^d mod p

	v[0] = g; v[1] = h; v[2] = p; v[3] = d;

	mpz_clear(q); mpz_clear(temp); gmp_randclear(state);
	return v;
}

std::vector<mpz_class> Encrypt(mpz_class m, mpz_class g, mpz_class h, mpz_class p)
{
	mpz_class c1, c2;
	mpz_t r; mpz_init(r);
	gmp_randstate_t state; gmp_randinit_mt(state); gmp_randseed_ui(state, (96889010407 * time(NULL) + 1013904223) % (unsigned long)pow(2, 31) - 1);
	mpz_t temp; mpz_init(temp); mpz_sub_ui(temp, p.get_mpz_t(), 2);

	mpz_urandomm(r, state, temp); mpz_add_ui(r, r, 1); //random integer r

	mpz_powm(c1.get_mpz_t(), g.get_mpz_t(), r, p.get_mpz_t()); //c1 = g^r mod p
	mpz_powm(c2.get_mpz_t(), h.get_mpz_t(), r, p.get_mpz_t()); mpz_mul(c2.get_mpz_t(), c2.get_mpz_t(), m.get_mpz_t()); mpz_mod(c2.get_mpz_t(), c2.get_mpz_t(), p.get_mpz_t()); //c2 = m*h^r mod p

	std::vector<mpz_class> v;
	v.emplace_back(c1); v.emplace_back(c2);

	mpz_clear(r); mpz_clear(temp); gmp_randclear(state);
	return v;
}

mpz_class Decrypt(mpz_class c1, mpz_class c2, mpz_class d, mpz_class p)
{
	mpz_class m;
	mpz_t s; mpz_init(s);

	mpz_powm(s, c1.get_mpz_t(), d.get_mpz_t(), p.get_mpz_t()); //s = c1^d mod p

	mpz_t inv; mpz_init(inv);
	mpz_invert(inv, s, p.get_mpz_t());

	mpz_mul(m.get_mpz_t(), c2.get_mpz_t(), inv); //m = c2*s^-1 mod p = c2*c1^-d mod p
	mpz_mod(m.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());

	mpz_clear(s); mpz_clear(inv);
	return m;
}