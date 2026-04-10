.PHONY: check build install install-clean test clean

check:
	cargo fmt
	cargo clippy --all-targets --all-features -- -D warnings

build: check
	cargo build --release

test: check
	cargo test

install: build
	@echo "Installing osprey to ~/.cargo/bin/"
	@cp target/release/osprey ~/.cargo/bin/osprey
	@echo "Installed: $$(osprey --version 2>/dev/null || echo 'osprey')"

install-clean: clean install

clean:
	cargo clean
