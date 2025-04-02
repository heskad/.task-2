#!/bin/bash

# Скрипт для установки всех conda-окружений из папки envs/
# Требует: Miniconda/Anaconda и менеджер mamba

# --- Проверка наличия conda ---
if ! command -v conda &> /dev/null; then
    echo "Ошибка: Conda не установлена. Сначала установите Miniconda/Anaconda."
    exit 1
fi

# --- Установка mamba (опционально) ---
if ! command -v mamba &> /dev/null; then
    echo "Установка mamba для ускорения..."
    conda install -n base -c conda-forge mamba -y
fi

# --- Директория с окружениями ---
ENVS_DIR="envs"
if [ ! -d "$ENVS_DIR" ]; then
    echo "Ошибка: Папка $ENVS_DIR не найдена. Создайте её и поместите туда .yaml-файлы."
    exit 1
fi

# --- Установка каждого окружения ---
for env_file in "$ENVS_DIR"/*.yaml; do
    if [ -f "$env_file" ]; then
        env_name=$(basename "$env_file" .yaml)
        echo "Устанавливаем окружение: $env_name из $env_file..."
        
        # Используем mamba, если доступна (иначе conda)
        if command -v mamba &> /dev/null; then
            mamba env create -f "$env_file" -n "$env_name" --force
        else
            conda env create -f "$env_file" -n "$env_name" --force
        fi

        if [ $? -eq 0 ]; then
            echo "Окружение $env_name успешно установлено."
        else
            echo "шибка при установке $env_name!"
        fi
    fi
done

echo "Готово! Все окружения установлены."