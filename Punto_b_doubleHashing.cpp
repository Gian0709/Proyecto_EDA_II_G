#ifndef PUNTO_B_DOUBLE_HASHING_H
#define PUNTO_B_DOUBLE_HASHING_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector> // Necesario para std::vector

/**
 * Estructura auxiliar para mapear un K-mer y su frecuencia acumulada
 * hacia estructuras externas como el Max-Heap.
 */
struct KmerFrequency {
    std::string kmer;
    int frequency;
};

/**
 * Estructura que representa cada celda individual (entrada) de la tabla hash.
 * Almacena la clave (el K-mer), el valor (frecuencia) y una bandera de borrado lógico.
 */
struct KmerEntry {
    std::string key;
    int value;
    bool isDeleted;

    // Constructor para inicializar una nueva entrada con valor inicial 1
    KmerEntry(const std::string& k, int v) {
        key = k;
        value = v;
        isDeleted = false;
    }
};

/**
 * Clase principal que gestiona la Tabla Hash con resolución de colisiones
 * mediante Doble Hashing.
 */
class DoubleHashTable {
private:
    KmerEntry** bins; // Arreglo dinámico de punteros a entradas (tabla de dispersión)
    int size;         // Tamaño total de la tabla (calculado típicamente como 4^n)
    int num_keys;     // Cantidad actual de elementos almacenados

    /**
     * Primera función hash basada en Horner's Method multiplicativo.
     */
    int hash1(const std::string& key) const {
        unsigned long hash = 0;
        for (size_t i = 0; i < key.length(); ++i) {
            hash = (hash * 31 + key[i]) % size;
        }
        return static_cast<int>(hash);
    }

    /**
     * Segunda función hash auxiliar utilizada para calcular el paso de desplazamiento
     * en caso de colisión (evita agrupamientos secundarios).
     */
    int hash2(const std::string& key) const {
        unsigned long hash = 0;
        for (size_t i = 0; i < key.length(); ++i) {
            hash = (hash * 37 + key[i]);
        }
        int step = 1 + (hash % (size - 1));
        return step;
    }

public:
    /**
     * Constructor: Inicializa la tabla dimensionada en función de la longitud del K-mer (4^n).
     */
    DoubleHashTable(int n) {
        size = static_cast<int>(std::pow(4, n));
        num_keys = 0;
        bins = new KmerEntry*[size];
        for (int i = 0; i < size; ++i) {
            bins[i] = nullptr;
        }
    }

    /**
     * Destructor: Libera toda la memoria dinámica reservada para las entradas y el arreglo de celdas.
     */
    ~DoubleHashTable() {
        for (int i = 0; i < size; ++i) {
            delete bins[i];
        }
        delete[] bins;
    }

    /**
     * Inserta un nuevo K-mer o incrementa su frecuencia si ya existe en la tabla.
     * Utiliza doble sondeo y maneja celdas previamente marcadas como borradas.
     */
    bool insertOrIncrement(const std::string& key) {
        int h1 = hash1(key);
        int h2 = hash2(key);
        int i = 0;

        int firstDeletedIndex = -1;

        while (i < size) {
            int index = (h1 + i * h2) % size;

            // 1. Si la celda está completamente vacía
            if (bins[index] == nullptr) {
                int targetIndex = (firstDeletedIndex != -1) ? firstDeletedIndex : index;

                if (firstDeletedIndex != -1) {
                    delete bins[targetIndex];
                }

                bins[targetIndex] = new KmerEntry(key, 1);
                num_keys++;
                return true;
            }

            // 2. Si la celda contiene exactamente la clave que se busca
            if (!bins[index]->isDeleted && bins[index]->key == key) {
                bins[index]->value++; // Incrementamos su frecuencia de aparición
                return true;
            }

            // 3. Si la celda está marcada como borrada, se registra la posición reutilizable
            if (bins[index]->isDeleted && firstDeletedIndex == -1) {
                firstDeletedIndex = index;
            }

            i++;
        }

        // Si se recorrió toda la tabla y se encontró una celda borrada reutilizable
        if (firstDeletedIndex != -1) {
            delete bins[firstDeletedIndex];
            bins[firstDeletedIndex] = new KmerEntry(key, 1);
            num_keys++;
            return true;
        }

        return false;
    }

    /**
     * Busca una clave (K-mer) en la tabla hash y retorna un puntero a su valor (frecuencia).
     * Retorna nullptr si el elemento no se encuentra.
     */
    int* lookup(const std::string& key) const {
        int h1 = hash1(key);
        int h2 = hash2(key);
        int i = 0;

        while (i < size) {
            int index = (h1 + i * h2) % size;
            if (bins[index] == nullptr) break;

            if (!bins[index]->isDeleted && bins[index]->key == key) {
                return &(bins[index]->value);
            }
            i++;
        }
        return nullptr;
    }

    // Métodos accesorios para consultar propiedades de la tabla
    int getSize() const { return size; }
    int getNumKeys() const { return num_keys; }


    // Muestra por consola todos los K-mers válidos y sus frecuencias almacenadas.

    void display() const {
        std::cout << "--- TABLA HASH 2 (DOUBLE HASHING) [Tamaño: " << size << "] ---\n";
        for (int i = 0; i < size; ++i) {
            if (bins[i] != nullptr && !bins[i]->isDeleted) {
                std::cout << "[" << i << "]: K-mer = " << bins[i]->key
                          << " | Frecuencia = " << bins[i]->value << "\n";
            }
        }
    }

    /**
     * Exporta todos los K-mers y frecuencias válidos hacia un vector estándar
     * para facilitar su procesamiento posterior en el Max-Heap.
     */
     // Codigo generado
    std::vector<KmerFrequency> obtenerKmersValidos() const {
        std::vector<KmerFrequency> lista;
        for (int i = 0; i < size; ++i) {
            if (bins[i] != nullptr && !bins[i]->isDeleted) {
                lista.push_back({bins[i]->key, bins[i]->value});
            }
        }
        return lista;
    }
};

#endif
