# Makefile для быстрых команд

.PHONY: install test clean run analyze visualize help

help:
	@echo "Доступные команды:"
	@echo "  make install     - Установка окружения"
	@echo "  make test        - Запуск тестов"
	@echo "  make run         - Запуск примера пайплайна"
	@echo "  make analyze     - Анализ результатов"
	@echo "  make visualize   - Визуализация"
	@echo "  make clean       - Очистка временных файлов"

install:
	conda env create -f environment.yml
	@echo "✓ Окружение создано. Активируй: conda activate protein_dynamics"

test:
	pytest tests/ -v --tb=short

test-fast:
	pytest tests/ -v -m "not slow"

run:
	python scripts/01_full_pipeline.py \
		--config config/pipeline_config.yaml \
		--output example_results

analyze:
	python scripts/02_analyze_results.py \
		--trajectory example_results/trajectory.dcd \
		--topology example_results/complex.pdb \
		--output example_results/analysis

visualize:
	python scripts/03_visualize.py \
		--pdb example_results/complex.pdb \
		--trajectory example_results/trajectory.dcd \
		--output example_results/viz \
		--create-gif

clean:
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	find . -type f -name "*.log" -delete
	rm -rf build/ dist/ *.egg-info/

clean-results:
	rm -rf results/ *_results/ output/ *_output/

format:
	black src/ scripts/ tests/
	isort src/ scripts/ tests/

lint:
	flake8 src/ scripts/ tests/
	mypy src/

docker-build:
	docker build -t protein-ligand-md .

docker-run:
	docker run -v $(PWD)/data:/data protein-ligand-md
