<!-- README.md is generated from README.Rmd. Please edit that file -->

# lkj

<!-- badges: start -->
<!-- Add badges here, e.g., CRAN, R-CMD-check, etc. -->
<!-- badges: end -->

The **`lkj`** package is an R tool designed for working with correlation matrices generated from the **Lewandowski–Kurowicka–Joe (LKJ)** distribution. It supports both **C-vine** and **onion** construction methods, providing functionalities for deterministic and random generation, optimizing tasks related to statistical analysis and simulation.

## Installation

To install the development version of the package from [GitHub](https://github.com), use:

```r
# Install from GitHub
install.packages("devtools")
devtools::install_github("alrobles/lkj")
```

> Currently not available on CRAN. If available in the future, instructions for installing from CRAN will be included.

---

## Main Features

- **Direct implementations in R and C++**:
  - Supports **C-vine** and **onion** methods.
- **Deterministic and Random Generation**:
  - Converts unconstrained reals to Cholesky factors.
  - Generates random correlation matrices.
- **Core Functions**:
  - `lkj()`: Deterministic mapping.
  - `rlkj()`: Random matrix generation.
- Mathematical documentation included alongside the source code.

---

## Quick Example

Here’s a quick example of how to use the package for deterministic mapping and random generation:

```r
library(lkj)

# Example 1: Deterministic mapping using C-vine and specific parameters
v <- c(-0.5, 0.3, 1.2, -0.8, 0.1, 0.6)
R <- lkj(v, d = 4, eta = 1, method = "cvine")
print(round(R, 4))

# Example 2: Random correlation matrix generation
R_random <- rlkj(10, d = 3, eta = 0.9)
print(R_random)
```

---

## Additional Documentation

For more information about this package (e.g., mathematical concepts involved or advanced examples), please explore the following vignettes from R:

```r
vignette("method-lkj")
vignette("roadmap")
```

You can also consult internal files such as:

- **Mathematical Methods**: `METHODS.md`.
- **Development Roadmap**: `roadmap.md`.

---

## License

The **`lkj`** package is licensed under **GPL (≥ 3)**. You can check the [LICENSE.md](LICENSE.md) file in this repository for more details.

---

For any inquiries, improvements, or issue reports, please check our [Issues section](https://github.com/alrobles/lkj/issues) in this repository.

Made with ❤️ by Angel Robles.

---

# lkj

<!-- badges: start -->
<!-- badges: end -->

El paquete **`lkj`** es una herramienta en R diseñada para trabajar con matrices de correlación generadas a partir de la distribución **Lewandowski–Kurowicka–Joe (LKJ)**. Ofrece soporte tanto para métodos de construcción **C-vine** como **onion**, con funcionalidades para generación determinística y aleatoria, optimizando tareas relacionadas con análisis estadístico y simulación.

## Instalación

Para instalar la versión en desarrollo del paquete desde [GitHub](https://github.com), usa:

```r
# Instalar desde GitHub
install.packages("devtools")
devtools::install_github("alrobles/lkj")
```

> Aún no disponible en CRAN. Si se encuentra disponible en el futuro, incluiremos las instrucciones pertinentes para instalar desde CRAN.

---

## Características principales

- **Implementaciones directas en R y C++**:
  - Métodos **C-vine** y **onion** para trabajar con matrices.
- **Soporte a generación determinística y aleatoria**:
  - Conversión de reals no restringidos a factores de Cholesky.
  - Generación de matrices aleatorias de correlación.
- **Funciones principales**:
  - `lkj()`: Mapeo determinista.
  - `rlkj()`: Generación aleatoria de matrices.
- Documentación matemática acompañada del código fuente.

---

## Ejemplo rápido

Aquí se muestra un ejemplo básico de uso del paquete para el mapeo determinista y la generación aleatoria:

```r
library(lkj)

# Ejemplo 1: Mapeo determinista usando C-vine y parámetros específicos
v <- c(-0.5, 0.3, 1.2, -0.8, 0.1, 0.6)
R <- lkj(v, d = 4, eta = 1, method = "cvine")
print(round(R, 4))

# Ejemplo 2: Generación aleatoria de matrices de correlación
R_random <- rlkj(10, d = 3, eta = 0.9)
print(R_random)
```

---

## Documentación adicional

Si necesitas más información sobre este paquete (por ejemplo, conceptos matemáticos involucrados o ejemplos más avanzados), te sugerimos explorar los siguientes manuales/vinetas desde R (presiona `Ctrl + Enter`):

```r
vignette("method-lkj")
vignette("roadmap")
```

También puedes consultar archivos internos como:

- **Métodos matemáticos**: `METHODS.md`.
- **Planes** de desarrollo: `roadmap.md`.

---

## Licencia

El paquete **`lkj`** está bajo la licencia **GPL (≥ 3)**. Puedes consultar el archivo [LICENSE.md](LICENSE.md) dentro de este repositorio para más detalles.

---

Con cualquier consulta, mejora o reporte de errores, por favor revisa [la sección Issues](https://github.com/alrobles/lkj/issues) en este repositorio.

Generado con ❤️ por Angel Robles.
