# Adleman–Pomerance–Rumely-Cohen-Lenstra primality test (APR-CL)

# Overview

Python implementation of APR-CL primality test based on [Cohen].

Primality test takes few seconds for 100 digit number, and 10 minutes for 300 digit number.



APR-CL法(Adleman–Pomerance–Rumely-Cohen-Lenstra)による素数判定プログラムです。Pythonで書かれています。

100桁なら数秒、300桁なら数分で素数判定ができます。



- 数学的な内容については、完全ではありませんが日本語で以下に解説したので、こちらも参照してください。詳細は論文を読んだ方が正確です。

https://wacchoz.hatenablog.com/entry/2018/12/16/155500 (in Japanese)



- Wikipedia

https://en.wikipedia.org/wiki/Adleman%E2%80%93Pomerance%E2%80%93Rumely_primality_test



# Reference

- [APR] Adleman, Leonard M.; Pomerance, Carl; Rumely, Robert S. (1983). "On distinguishing prime numbers from composite numbers". Annals of Mathematics. 117 (1): 173–206.
- [Coh-Len1] Cohen, Henri; Lenstra, Hendrik W., Jr. (1984). "Primality testing and Jacobi sums". Mathematics of Computation. 42 (165): 297–330.
- [Coh-Len2] H. Cohen and A. K. Lenstra, Implementation of a new primality test. Math. Comp. 48 (1987), 103-121.
- [Cohen] H. Cohen, "A course in computational algebraic number theory". Springer-Verlag, Berlin, 1993