.PHONY: test

install:
	julia -e 'using Pkg; Pkg.develop(path=".");'

uninstall:
	julia -e 'using Pkg; Pkg.rm("MatricesForHomalg");'

test:
	julia -e 'using Pkg; Pkg.test("MatricesForHomalg");'
