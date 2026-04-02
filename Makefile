.PHONY: check build install test

check:
	cargo fmt
	cargo clippy --all-targets --all-features -- -D warnings

build: check
	cargo build --release

test: check
	cargo test

install: build
	cargo install --path crates/osprey
