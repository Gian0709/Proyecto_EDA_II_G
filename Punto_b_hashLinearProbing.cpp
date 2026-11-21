#ifndef PUNTO_B_LINEAR_PROBING_H
#define PUNTO_B_LINEAR_PROBING_H

#include <iostream>
#include <string>
#include <vector>

/**
 * Estructura que representa cada celda individual (entrada) de la tabla hash.
 * Almacena la clave (ID), el valor (Secuencia de ADN) y una bandera de borrado lógico.
 */
struct Entry {
    std::string key;
    std::string value;
    bool isDeleted;

    // Constructor para inicializar una nueva entrada con borrado lógico en falso
    Entry(const std::string& k, const std::string& v) {
        key = k;
        value = v;
        isDeleted = false;
    }
};

/**
 * Clase principal que gestiona la Tabla Hash con resolución de colisiones
 * mediante Sondeo Lineal y redimensionamiento automático.
 */
class HashTable {
private:
    Entry** bins;           // Arreglo dinámico de punteros a entradas (tabla de dispersión)
    int size;               // Tamaño actual de la tabla
    int num_keys;           // Cantidad actual de elementos almacenados
    static const int DEFAULT_SIZE = 10; // Tamaño inicial por defecto

    /**
     * Función auxiliar para verificar si un número entero es primo.
     */
    bool isPrime(int n) const {
        if (n <= 1) return false;
        if (n <= 3) return true;
        if (n % 2 == 0 || n % 3 == 0) return false;
        for (int i = 5; i * i <= n; i += 6) {
            if (n % i == 0 || n % (i + 2) == 0) return false;
        }
        return true;
    }

    /**
     * Función auxiliar para encontrar el siguiente número primo a partir de 'n'
     * (útil para dimensionar la tabla hash con tamaños primos).
     */
    int nextPrime(int n) const {
        if (n <= 1) return 2;
        int prime = n;
        while (!isPrime(prime)) {
            prime++;
        }
        return prime;
    }

    /**
     * Función hash basada en Horner's Method multiplicativo para transformar claves en índices.
     */
    int hashFunction(const std::string& key) const {
        unsigned long hash = 0;
        for (size_t i = 0; i < key.length(); ++i) {
            hash = (hash * 31 + key[i]) % size;
        }
        return static_cast<int>(hash);
    }

    /**
     * Método privado para redimensionar y reestructurar la tabla cuando
     * el factor de carga supera o iguala el límite (Load Factor >= 1.0).
     */
    void resize() {
        int oldSize = size;
        Entry** oldBins = bins;

        // Nuevo tamaño: siguiente número primo que duplica el tamaño actual
        size = nextPrime(oldSize * 2);
        bins = new Entry*[size];
        for (int i = 0; i < size; ++i) {
            bins[i] = nullptr;
        }

        num_keys = 0; // Se reinicia el contador para reinsertar limpiamente

        // Reinsertar los elementos válidos de la tabla antigua a la nueva
        for (int i = 0; i < oldSize; ++i) {
            if (oldBins[i] != nullptr && !oldBins[i]->isDeleted) {
                insertWithoutResize(oldBins[i]->key, oldBins[i]->value);
                delete oldBins[i];
            } else if (oldBins[i] != nullptr) {
                delete oldBins[i];
            }
        }
        delete[] oldBins;
    }

    /**
     * Inserción interna auxiliar sin evaluar redimensionamiento
     * (empleada exclusivamente durante el proceso de resize).
     */
    bool insertWithoutResize(const std::string& key, const std::string& value) {
        int index = hashFunction(key);
        int startIndex = index;

        // Sondeo lineal para buscar espacio disponible o clave existente
        while (bins[index] != nullptr && !bins[index]->isDeleted && bins[index]->key != key) {
            index = (index + 1) % size;
            if (index == startIndex) return false;
        }

        if (bins[index] != nullptr && !bins[index]->isDeleted && bins[index]->key == key) {
            bins[index]->value = value;
            return true;
        }

        if (bins[index] != nullptr && bins[index]->isDeleted) {
            delete bins[index];
            bins[index] = nullptr;
        }

        bins[index] = new Entry(key, value);
        num_keys++;
        return true;
    }

public:
    /**
     * Constructor: Inicializa la tabla hash con un tamaño predeterminado o personalizado.
     */
    HashTable(int s = DEFAULT_SIZE) : size(s), num_keys(0) {
        bins = new Entry*[size];
        for (int i = 0; i < size; ++i)
            bins[i] = nullptr;
    }

    /**
     * Destructor: Libera toda la memoria dinámica asignada a las entradas y al arreglo de celdas.
     */
    ~HashTable() {
        for (int i = 0; i < size; ++i)
            delete bins[i];
        delete[] bins;
    }

    /**
     * Método público de inserción unificado y seguro.
     * Evalúa el factor de carga y redimensiona automáticamente de ser necesario.
     */
    bool insert(const std::string& key, const std::string& value) {
        // 1. Verificar factor de carga: si es >= 1.0, ampliamos la tabla primero (corrección aplicada)
        double loadFactor = static_cast<double>(num_keys) / size;
        if (loadFactor >= 1.0) {
            resize();
        }

        // 2. Aplicar función hash y sondeo lineal con el tamaño actualizado
        int index = hashFunction(key);
        int startIndex = index;

        while (bins[index] != nullptr && !bins[index]->isDeleted && bins[index]->key != key) {
            index = (index + 1) % size;
            if (index == startIndex) return false; // Tabla llena
        }

        // Si la clave ya se encuentra registrada, actualizamos su valor asociado
        if (bins[index] != nullptr && !bins[index]->isDeleted && bins[index]->key == key) {
            bins[index]->value = value;
            return true;
        }

        // Si encontramos una celda nula o marcada previamente como borrada (dummy)
        if (bins[index] == nullptr || bins[index]->isDeleted) {
            if (bins[index] != nullptr && bins[index]->isDeleted) {
                delete bins[index];
            }
            bins[index] = new Entry(key, value);
            num_keys++;
            return true;
        }

        return false;
    }

    /**
     * Busca una clave (ID) en la tabla hash mediante sondeo lineal.
     * Retorna un puntero al valor asociado o nullptr si no se encuentra.
     */
    std::string* lookup(const std::string& key) const {
        int index = hashFunction(key);
        int count = 0;
        while (bins[index] != nullptr && count < size) {
            if (!bins[index]->isDeleted && bins[index]->key == key) {
                return &(bins[index]->value);
            }
            index = (index + 1) % size;
            count++;
        }
        return nullptr;
    }

    /**
     * Realiza un borrado lógico de una clave marcándola con la bandera 'isDeleted'.
     */
    bool remove(const std::string& key) {
        int index = hashFunction(key);
        int count = 0;
        while (bins[index] != nullptr && count < size) {
            if (!bins[index]->isDeleted && bins[index]->key == key) {
                bins[index]->isDeleted = true;
                num_keys--;
                return true;
            }
            index = (index + 1) % size;
            count++;
        }
        return false;
    }

    /**
     * Imprime el estado actual de la tabla hash (tamaño, elementos y contenido de cada celda).
     */
    void display() const {
        std::cout << "Tamaño actual de tabla: " << size << " | Elementos: " << num_keys << "\n";
        for (int i = 0; i < size; ++i) {
            std::cout << "[" << i << "]: ";
            if (bins[i] == nullptr) {
                std::cout << "null";
            } else if (bins[i]->isDeleted) {
                std::cout << "deleted";
            } else {
                std::cout << "(" << bins[i]->key << ", " << bins[i]->value << ")";
            }
            std::cout << "\n";
        }
    }
};

#endif
