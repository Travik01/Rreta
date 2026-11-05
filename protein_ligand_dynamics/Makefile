.PHONY: help install test test-fast run analyze visualize clean clean-results docker-build docker-run

help:
	@echo "Доступные команды:"
	@echo "  make install       - создать Conda-окружение"
	@echo "  make test          - запустить тесты"
	@echo "  make test-fast     - запустить быстрые тесты"
	@echo "  make run           - запустить пайплайн"
	@echo "  make analyze       - проанализировать траекторию"
	@echo "  make visualize     - построить визуализации"
	@echo "  make clean         - удалить кэш и сборки"
	@echo "  make clean-results - удалить каталоги результатов"
	@echo "  make docker-build  - собрать образ"
	@echo "  make docker-run    - запустить контейнер"

install:
	conda env create -f environment.yml

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

docker-build:
	docker build -t protein-ligand-md .

docker-run:
	docker run -v $(PWD)/results:/app/results protein-ligand-md
