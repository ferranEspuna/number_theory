def generate_primes(limit=float('inf')):
    primes = [2]
    n = 3
    while True:
        if n > limit:
            break
        for p in primes:
            if n % p == 0:
                break
            if p * p >= n:
                primes.append(n)
                yield n
                break
        n += 1