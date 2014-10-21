import random

cache = {}

def count_partitions(n, limit):
    if n == 0:
        return 1
    if (n, limit) in cache:
        return cache[n, limit]
    x = cache[n, limit] = sum(count_partitions(n-k, k) for k in range(1, min(limit, n) + 1))
    return x

def random_partition(n):
    a = []
    limit = n
    total = count_partitions(n, limit)
    which = random.randrange(total)
    while n:
        for k in range(1, min(limit, n) + 1):
            count = count_partitions(n-k, k)
            if which < count:
                break
            which -= count
        a.append(k)
        limit = k
        n -= k
    return a

