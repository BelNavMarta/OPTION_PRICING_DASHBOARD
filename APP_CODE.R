############################################################################## 
# ÍNDICE DEL CÓDIGO:
#   1. Librerías
#   2. Funciones matemáticas (Binomial + Black-Scholes + Griegas + Monte Carlo)
#   3. UI (Interfaz de usuario)
#   4. Server (Lógica reactiva)
#   5. Lanzamiento de la app
############################################################################### 


############################################################################### 
# 1. LIBRERÍAS
############################################################################### 

#install.packages("shiny")
#install.packages("plotly")

library(shiny)      
library(plotly)      # Gráficos interactivos

# Si x no existe, es NULL o está vacío, usa y. Si x sí existe, usa x:
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

############################################################################### 
# 2. FUNCIONES MATEMÁTICAS
############################################################################### 

# -----------------------------------------------------------------------------
# 2.1 MODELO BINOMIAL 
# -----------------------------------------------------------------------------

# (Wilmott, Cap. 3, p. 60-86):

# El modelo binomial discretiza el movimiento del precio del subyacente S.
# En cada paso temporal δt, el precio puede:
#   - Subir:  S → u·S   con probabilidad risk-neutral p
#   - Bajar:  S → d·S   con probabilidad risk-neutral 1-p

# PARAMETRIZACIÓN:
#   δt = T / N  → tamaño del paso temporal
#   u  = exp(σ · √δt) → factor de subida
#   d  = 1/u = exp(-σ · √δt)  → factor de bajada
#   p  = (exp((r-q)·δt) - d) / (u - d) → probabilidad riesgo-neutral

# ÁRBOL DE PRECIOS (nodo n, nivel j desde abajo):
#   S[n,j] = S₀ · u^j · d^(n-j),  con j = 0, 1, ..., n

# VALORACIÓN (backward induction):
#   En T (vencimiento): V[N,j] = Payoff(S[N,j])
#   Hacia atrás: V[n,j] = e^(-r·δt) · [p·V[n+1,j+1] + (1-p)·V[n+1,j]]

# EJERCICIO ANTICIPADO (Americanas, Wilmott p. 88-89):
#   V[n,j] = max(Payoff(S[n,j]),  e^(-r·δt) · [p·V[n+1,j+1] + (1-p)·V[n+1,j]])
#   En cada nodo comprobamos si ejercer ya es mejor que esperar


binomial_model <- function(S, K, r, q, sigma, T, N, tipo_opcion, estilo_ejercicio) {
  # Parámetros de entrada:
  #   S: precio actual del subyacente (S₀)
  #   K: precio de ejercicio
  #   r: tasa libre de riesgo continua anualizada
  #   q: tasa continua anualizada de dividendos del subyacente
  #   sigma: volatilidad anualizada del subyacente
  #   T: tiempo hasta vencimiento en años
  #   N: número de pasos del árbol
  #   tipo_opcion: call o put
  #   estilo_ejercicio: european o american
  
  
  # Cálculo de los parámetros del árbol:
  
  dt <- T / N   # δt: tamaño de cada paso temporal 
  
  u  <- exp(sigma * sqrt(dt))   # factor de subida 
  d  <- 1 / u                   # factor de bajada; d = 1/u garantiza árbol recombinante: u·d = 1, por lo que subir y luego bajar = bajar y luego subir → los nodos centrales coinciden
  
  # Probabilidad riesgo-neutral (Wilmott, p. 72-74):
  # Se obtiene de: exp((r-q)·δt) = p·u + (1-p)·d → despejando p:
  p <- (exp((r - q) * dt) - d) / (u - d)    
  
  # Factor de descuento por paso (Wilmott, p. 80):
  # Descontar un paso temporal al tipo libre de riesgo
  descuento <- exp(-r * dt)
  
  
  # Construimos el árbol de precios del subyacente:
  # S[n,j] = S₀ · u^j · d^(n-j)
  # Guardamos todos los nodos en una matriz para la visualización del árbol
  # Filas (n) = paso temporal (0..N), Columnas (j) = nivel del nodo (0..N) 
  
  S_tree <- matrix(NA, nrow = N + 1, ncol = N + 1)
  for (n in 0:N) {
    for (j in 0:n) {
      # En el paso n, el nodo j (desde abajo) tiene precio:
      S_tree[n + 1, j + 1] <- S * (u^j) * (d^(n - j))
    }
  }
  
  
  # Calculamos payoffs en el vencimiento (n = N):
  #   Call: max(S - K, 0)
  #   Put:  max(K - S, 0)
  
  V_tree <- matrix(NA, nrow = N + 1, ncol = N + 1)
  intrinsic_tree <- matrix(0, nrow = N + 1, ncol = N + 1)
  continuation_tree <- matrix(NA, nrow = N + 1, ncol = N + 1)
  exercise_tree <- matrix(FALSE, nrow = N + 1, ncol = N + 1)
  
  for (j in 0:N) {
    S_final <- S_tree[N + 1, j + 1]
    if (tipo_opcion == "call") {
      V_tree[N + 1, j + 1] <- max(S_final - K, 0)
    } else {  # put
      V_tree[N + 1, j + 1] <- max(K - S_final, 0)
    }
    intrinsic_tree[N + 1, j + 1] <- V_tree[N + 1, j + 1]
  }
  
  
  # Backward induction (recorremos el árbol hacia atrás):
  # Partimos del vencimiento y nos movemos hacia el presente.
  # Wilmott (p. 80-86): 
  
  for (n in (N - 1):0) {          # recorremos de N-1 hasta 0
    for (j in 0:n) {               # en el paso n hay n+1 nodos
      
      # Valor de continuar (esperar a ejercer): Fórmula de pricing binomial
      V_continuar <- descuento * (p * V_tree[n + 2, j + 2] +
                                    (1 - p) * V_tree[n + 2, j + 1])
      
      # Valor de ejercer inmediatamente (solo importa para Americanas):
      S_nodo <- S_tree[n + 1, j + 1]
      if (tipo_opcion == "call") {
        V_ejercer <- max(S_nodo - K, 0)
      } else {
        V_ejercer <- max(K - S_nodo, 0)
      }
      
      continuation_tree[n + 1, j + 1] <- V_continuar
      intrinsic_tree[n + 1, j + 1] <- V_ejercer
      
      if (estilo_ejercicio == "american") {
        # Para americanas: elegimos el máximo entre ejercer y esperar
        # Wilmott (p. 88-89): "we must ensure there are no arbitrage opportunities at any of the nodes"
        V_tree[n + 1, j + 1] <- max(V_continuar, V_ejercer)
        # Marcamos ejercicio anticipado solo si ejercer ahora es estrictamente mejor que continuar.
        # La tolerancia evita marcar nodos por pequeñas diferencias numéricas.
        exercise_tree[n + 1, j + 1] <- (V_ejercer > V_continuar + 1e-10)  
      } else {
        # Para europeas: solo podemos esperar (no hay ejercicio anticipado)
        V_tree[n + 1, j + 1] <- V_continuar
        exercise_tree[n + 1, j + 1] <- FALSE
      }
    }
  }
  
  # El precio de la opción hoy es el valor en el nodo raíz (n=0, j=0)
  precio_opcion <- V_tree[1, 1]
  
  # Devolvemos también los árboles completos para visualización y griegas
  return(list(
    precio   = precio_opcion,
    S_tree   = S_tree,
    V_tree   = V_tree,
    intrinsic_tree = intrinsic_tree,
    continuation_tree = continuation_tree,
    exercise_tree = exercise_tree,
    u        = u,
    d        = d,
    p        = p,
    dt       = dt,
    N        = N
  ))
}


# MODELO BINOMIAL PARA BARRERAS 

# (Wilmott, Cap. 13)

# En cada nodo del árbol se puede llegar por caminos vivos o caminos muertos.
#   - Para Knock-Out: "vivo" = el precio nunca ha cruzado la barrera.
#   - Para Knock-In:  "vivo" = el precio sí ha cruzado la barrera (activado).

# Si no hay ningún camino vivo que llegue a un nodo, el nodo es inactivo y
# V = 0 siempre (para knock-out) o en el caso knock-in la valoración continúa
# porque caminos futuros pueden activar la opción.

# ALGORITMO KNOCK-OUT:
#   1. alive[n,j] = (S[n,j] no cruza barrera) y (algún padre es alive)
#   2. Nodo terminal: V[N,j] = payoff si alive[N,j], 0 si no.
#   3. Backward: V[n,j] = 0 si not alive; fórmula estándar si alive.
#      (Americana: max(cont, ejercer) si alive)

# ALGORITMO KNOCK-IN:
#   1. activated[n,j] = (S[n,j] cruza barrera) o (algún padre activado)
#   2. Terminal: V[N,j] = payoff si activated[N,j], 0 si no.
#   3. Backward estándar para todos los nodos (incluso no activados pueden
#      tener valor positivo si futuros caminos activan la barrera).
#      (American: max(cont, ejercer) solo si activated[n,j])

# INDEXACIÓN DEL ÁRBOL:
#   Nodo (step=n, level=j) → S_tree[n+1, j+1]
#   Padres de (n, j):
#     - Por movimiento UP desde (n-1, j-1): existe si j >= 1
#     - Por movimiento DOWN desde (n-1, j): existe si j <= n-1
# -----------------------------------------------------------------------------

binomial_barrier_model <- function(S, K, r, q, sigma, T, N, tipo_opcion,
                                   estilo_ejercicio, tipo_barrera, B) {
  dt <- T / N
  u  <- exp(sigma * sqrt(dt))
  d  <- 1 / u
  p  <- (exp((r - q) * dt) - d) / (u - d)
  desc <- exp(-r * dt)
  
  # Árbol de precios (idéntico al vanilla)
  S_tree <- matrix(NA, nrow = N + 1, ncol = N + 1)
  for (n in 0:N) {
    for (j in 0:n) {
      S_tree[n + 1, j + 1] <- S * (u^j) * (d^(n - j))
    }
  }
  
  es_ko <- tipo_barrera %in% c("down-and-out", "up-and-out")
  # TRUE si la barrera es inferior: se toca cuando S_nodo <= B
  # FALSE si la barrera es superior: se toca cuando S_nodo >= B
  es_down <- tipo_barrera %in% c("down-and-out", "down-and-in") 
  
  barrera_tocada <- function(S_nodo) { 
    if (es_down) S_nodo <= B else S_nodo >= B
  }
  
  payoff <- function(S_nodo) {
    if (tipo_opcion == "call") max(S_nodo - K, 0) else max(K - S_nodo, 0)
  }
  
  intrinsic_tree <- matrix(0, nrow = N + 1, ncol = N + 1)
  for (n in 0:N) {
    for (j in 0:n) {
      intrinsic_tree[n + 1, j + 1] <- payoff(S_tree[n + 1, j + 1])
    }
  }
  
  exercise_tree <- matrix(FALSE, nrow = N + 1, ncol = N + 1)
  continuation_tree <- matrix(NA, nrow = N + 1, ncol = N + 1)
  
  if (es_ko) {
    # Knock-out: se valora el estado vivo. Si la barrera se toca, el valor pasa a cero.
    V_live <- matrix(0, nrow = N + 1, ncol = N + 1)
    alive_possible <- matrix(FALSE, nrow = N + 1, ncol = N + 1)
    alive_possible[1, 1] <- !barrera_tocada(S_tree[1, 1])
    
    for (n in 1:N) {
      for (j in 0:n) {
        if (barrera_tocada(S_tree[n + 1, j + 1])) {
          alive_possible[n + 1, j + 1] <- FALSE
        } else {
          parent_up   <- (j >= 1)   && alive_possible[n, j]
          parent_down <- (j <= n-1) && alive_possible[n, j + 1]
          alive_possible[n + 1, j + 1] <- parent_up || parent_down
        }
      }
    }
    
    for (j in 0:N) {
      V_live[N + 1, j + 1] <- if (alive_possible[N + 1, j + 1]) {
        payoff(S_tree[N + 1, j + 1])
      } else {
        0
      }
    }
    
    for (n in (N - 1):0) {
      for (j in 0:n) {
        
        if (!alive_possible[n + 1, j + 1]) {
          
          # Si no existe ningún camino vivo que llegue al nodo,
          # la opción knock-out ya está eliminada.
          V_live[n + 1, j + 1] <- 0
          continuation_tree[n + 1, j + 1] <- 0
          
        } else {
          
          # Si el hijo de subida no es alcanzable por un camino vivo, su valor es 0.
          up_value <- if (alive_possible[n + 2, j + 2]) {
            V_live[n + 2, j + 2]
          } else {
            0
          }
          
          # Si el hijo de bajada no es alcanzable por un camino vivo, su valor es 0.
          down_value <- if (alive_possible[n + 2, j + 1]) {
            V_live[n + 2, j + 1]
          } else {
            0
          }
          
          V_cont <- desc * (p * up_value + (1 - p) * down_value)
          continuation_tree[n + 1, j + 1] <- V_cont
          
          V_ex <- payoff(S_tree[n + 1, j + 1])
          
          if (estilo_ejercicio == "american" && V_ex > V_cont + 1e-10) {
            V_live[n + 1, j + 1] <- V_ex
            exercise_tree[n + 1, j + 1] <- V_ex > 0
          } else {
            V_live[n + 1, j + 1] <- V_cont
          }
        }
      }
    }
    
    V_tree <- V_live
    state_tree <- alive_possible
    V_inactive <- matrix(0, nrow = N + 1, ncol = N + 1)
    
  } else {
    # Knock-in: se necesitan dos estados por nodo.
    # V_not: valor si la opción todavía no ha sido activada.
    # V_act: valor si la opción ya ha sido activada.
    V_not <- matrix(0, nrow = N + 1, ncol = N + 1)
    V_act <- matrix(0, nrow = N + 1, ncol = N + 1)
    act_possible <- matrix(FALSE, nrow = N + 1, ncol = N + 1)
    not_possible <- matrix(FALSE, nrow = N + 1, ncol = N + 1)
    
    act_possible[1, 1] <- barrera_tocada(S_tree[1, 1])
    not_possible[1, 1] <- !barrera_tocada(S_tree[1, 1])
    
    for (n in 1:N) {
      for (j in 0:n) {
        hit <- barrera_tocada(S_tree[n + 1, j + 1])
        parent_act <- ((j >= 1) && act_possible[n, j]) || ((j <= n-1) && act_possible[n, j + 1])
        parent_not <- ((j >= 1) && not_possible[n, j]) || ((j <= n-1) && not_possible[n, j + 1])
        act_possible[n + 1, j + 1] <- hit || parent_act || (hit && parent_not)
        not_possible[n + 1, j + 1] <- (!hit) && parent_not
      }
    }
    
    for (j in 0:N) {
      V_act[N + 1, j + 1] <- payoff(S_tree[N + 1, j + 1])
      V_not[N + 1, j + 1] <- if (barrera_tocada(S_tree[N + 1, j + 1])) payoff(S_tree[N + 1, j + 1]) else 0
    }
    
    for (n in (N - 1):0) {
      for (j in 0:n) {
        # Estado activado: se comporta como una vanilla desde ese nodo.
        V_cont_act <- desc * (p * V_act[n + 2, j + 2] + (1 - p) * V_act[n + 2, j + 1])
        continuation_tree[n + 1, j + 1] <- V_cont_act
        V_ex <- payoff(S_tree[n + 1, j + 1])
        if (estilo_ejercicio == "american" && V_ex > V_cont_act + 1e-10) {
          V_act[n + 1, j + 1] <- V_ex
          exercise_tree[n + 1, j + 1] <- V_ex > 0
        } else {
          V_act[n + 1, j + 1] <- V_cont_act
        }
        
        # Estado no activado: al pasar a un hijo, puede activarse si el hijo toca la barrera.
        up_value <- if (barrera_tocada(S_tree[n + 2, j + 2])) V_act[n + 2, j + 2] else V_not[n + 2, j + 2]
        down_value <- if (barrera_tocada(S_tree[n + 2, j + 1])) V_act[n + 2, j + 1] else V_not[n + 2, j + 1]
        V_not[n + 1, j + 1] <- desc * (p * up_value + (1 - p) * down_value)
        if (!act_possible[n + 1, j + 1]) {
          continuation_tree[n + 1, j + 1] <- V_not[n + 1, j + 1]
        }
        
        # Si el nodo actual ya toca la barrera, el estado correcto es activado.
        if (barrera_tocada(S_tree[n + 1, j + 1])) {
          V_not[n + 1, j + 1] <- V_act[n + 1, j + 1]
        }
      }
    }
    
    V_tree <- if (barrera_tocada(S_tree[1, 1])) V_act else V_not
    state_tree <- act_possible
    V_inactive <- V_not
    V_active <- V_act
  }
  
  precio <- if (es_ko) V_tree[1, 1] else {
    if (barrera_tocada(S_tree[1, 1])) V_active[1, 1] else V_inactive[1, 1]
  }
  
  list(
    precio = precio,
    S_tree = S_tree,
    V_tree = if (es_ko) V_tree else ifelse(state_tree, V_active, V_inactive),
    V_live = if (es_ko) V_tree else NULL,
    V_active = if (!es_ko) V_active else NULL,
    V_not_active = if (!es_ko) V_inactive else NULL,
    intrinsic_tree = intrinsic_tree,
    continuation_tree = continuation_tree,
    exercise_tree = exercise_tree,
    alive_matrix = state_tree,
    es_ko = es_ko,
    u = u, d = d, p = p, dt = dt, N = N, B = B, tipo_barrera = tipo_barrera
  )
}


# -----------------------------------------------------------------------------
# 2.2 GRIEGAS BINOMIALES (aproximaciones numéricas)
# -----------------------------------------------------------------------------

# (Wilmott, Cap. 3, pp. 86-88):

# Las griegas miden la sensibilidad del precio de la opción a cambios en los
# parámetros. En el árbol binomial se calculan por diferencias finitas:

# DELTA (Δ): sensibilidad al precio del subyacente
# Ratio de cobertura: cuántas acciones necesitamos para cubrir la opción.
# Suponemos que hemos vendido una call. Si el precio sube, la call vale más y perdemos. 
# Para cubrirnos, compramos Δ acciones: si S sube 1€, la acción sube 1€ pero la call sube Δ€, 
# y si tenemos Δ acciones, nuestra ganancia en acciones cancela exactamente nuestra pérdida en la call.
# Wilmott: "delta hedging" - construir una cartera sin riesgo instantáneo, lo cual obliga al precio de la opción 
# a satisfacer la ecuación de Black-Scholes (argumento de no arbitraje). - el punto de partida de todo el modelo 

# GAMMA (Γ): curvatura; sensibilidad del delta al precio 
#   El delta nos dice cuánto vale la opción si S se mueve un poco. Pero Delta en sí cambia cuando S cambia (no es constante). 
#   Gamma mide esa variación de Delta. 
#   Una opicón con gamma alto necesita rehedgearse con mucha frecuencia porque su delta cambia rápido con el precio.

# THETA (Θ): sensibilidad al tiempo (decaimiento temporal)

# Delta, Gamma y Theta son derivadas respecto a variables de estado del árbol (S y t) — variables que ya están "dentro" del árbol.
# Vega y Rho son derivadas respecto a parámetros del modelo (σ y r) — valores que fijamos antes de construir el árbol.

# Por eso no podemos calcularlos mirando nodos vecinos. 
# En su lugar, necesitamos recalcular el árbol entero dos veces (perturbando σ o r), y tomar la diferencia.

# VEGA (ν): sensibilidad a la volatilidad 
#   No es derivada respecto a una variable del árbol, sino a un parámetro.
#   Se calcula con dos valoraciones binomiales con σ±ε:

# RHO (ρ): sensibilidad al tipo de interés 

# -----------------------------------------------------------------------------
# Las griegas se calculan sobre el árbol ya construido.
# Wilmott distingue dos tipos:

#   Tipo 1 — derivadas respecto a variables de estado (S y t):
#     Delta, Gamma, Theta: se obtienen directamente de los nodos del árbol,
#       sin necesidad de reconstruirlo.

#   Tipo 2 — derivadas respecto a parámetros del modelo (σ y r):
#     Vega, Rho: requieren reconstruir el árbol dos veces con el parámetro
#       perturbado, porque σ y r determinan la geometría del árbol (u, d, p).
# -----------------------------------------------------------------------------

griegas_binomiales <- function(resultado, S, K, r, q, sigma, T, N,
                               tipo_opcion, estilo_ejercicio) {
  
  S_tree <- resultado$S_tree
  V_tree <- resultado$V_tree
  dt     <- resultado$dt
  u      <- resultado$u      # factor de subida
  d      <- resultado$d      # factor de bajada
  
  # (Fig. 3.34 de Wilmott):
  V_plus  <- V_tree[2, 2]   # V⁺ — valor si precio subió (n=1, j=1)
  V_minus <- V_tree[2, 1]   # V⁻ — valor si precio bajó  (n=1, j=0)
  
  # ---------------------------------------------------------------------------
  # DELTA  (Wilmott, p. 86, Fig. 3.33)
  
  #          V⁺ − V⁻
  #   Δ  =  ─────────
  #          (u − d)·S
  
  # Ratio de cobertura: número de acciones del subyacente necesarias para
  # cubrir instantáneamente la posición en la opción (delta hedging).
  # En el límite δt → 0 converge a ∂V/∂S.
  # ---------------------------------------------------------------------------
  delta <- (V_plus - V_minus) / ((u - d) * S)
  
  
  # ---------------------------------------------------------------------------
  # GAMMA  (Wilmott, Fig. 3.34)
  
  # Γ = ∂²V/∂S² ≈ (Δ⁺ − Δ⁻) / distancia entre los puntos G
  
  # Se calculan dos deltas parciales usando los tres nodos del paso n=2
  # (marcados "G" en la figura):
  #   Δ⁺ = delta entre nodos (n=2,j=2) y (n=2,j=1)
  #   Δ⁻ = delta entre nodos (n=2,j=1) y (n=2,j=0)
  
  # La distancia entre los puntos medios donde se evalúan Δ⁺ y Δ⁻ es:
  #   ½·(S[2,2] − S[2,0]) = ½·(S·u² − S·d²)
  # ---------------------------------------------------------------------------
  if (N >= 2) {
    delta_up   <- (V_tree[3, 3] - V_tree[3, 2]) /
      (S_tree[3, 3] - S_tree[3, 2])   # Δ⁺ (nodos superiores de G)
    
    delta_down <- (V_tree[3, 2] - V_tree[3, 1]) /
      (S_tree[3, 2] - S_tree[3, 1])   # Δ⁻ (nodos inferiores de G)
    
    gamma <- (delta_up - delta_down) /
      (0.5 * (S_tree[3, 3] - S_tree[3, 1]))
  } else {
    gamma <- NA
  }
  
  
  # ---------------------------------------------------------------------------
  # THETA  
  
  # En CRR, u·d = 1 → el nodo central n=2, j=1 tiene precio exactamente S.
  # Comparamos V al mismo precio pero 2δt más tarde:
  
  #   Θ = (V[n=2, j=1] − V[n=0, j=0]) / (2·δt)
  
  # Esto mide el decaimiento temporal con el subyacente constante.
  # La fórmula ½(V⁺+V⁻) de Wilmott asume p≈½ (r≈0) y falla para r>0.
  # ---------------------------------------------------------------------------
  if (N >= 2) {
    theta        <- (V_tree[3, 2] - V_tree[1, 1]) / (2 * dt)
    theta_diario <- theta / 365 # convertido a €/día
  } else {
    theta        <- NA
    theta_diario <- NA
  }
  
  
  # ---------------------------------------------------------------------------
  # VEGA  (Wilmott, p. 88)
  
  #   ν  ≈  (V₊ − V₋) / (2ε)
  
  # Donde V₊ = V(σ+ε) y V₋ = V(σ−ε) son precios obtenidos reconstruyendo
  # el árbol completo con la volatilidad perturbada.
  # Nota de Wilmott: el precio base mejorado sería V ≈ ½(V₊ + V₋).
  # ---------------------------------------------------------------------------
  eps_vega     <- 0.001    # ε = 0.1% de volatilidad 
  V_sigma_up   <- binomial_model(S, K, r, q, sigma + eps_vega, T, N,
                                 tipo_opcion, estilo_ejercicio)$precio   # V₊
  V_sigma_down <- binomial_model(S, K, r, q, sigma - eps_vega, T, N,
                                 tipo_opcion, estilo_ejercicio)$precio   # V₋
  
  vega     <- (V_sigma_up - V_sigma_down) / (2 * eps_vega)
  vega_pct <- vega / 100   # por 1% de cambio en σ
  
  
  # ---------------------------------------------------------------------------
  # RHO
  
  #   ρ  ≈  (V₊ − V₋) / (2ε)
  
  # Misma lógica que vega: se reconstruye el árbol perturbando r en ±ε.
  # El tipo de interés afecta tanto a la probabilidad riesgo-neutral p
  # como al factor de descuento e^{−rδt}, por eso requiere re-valoración.
  # ---------------------------------------------------------------------------
  eps_rho    <- 0.0001    # ε = 0.01%
  V_r_up     <- binomial_model(S, K, r + eps_rho, q, sigma, T, N,
                               tipo_opcion, estilo_ejercicio)$precio   # V₊
  V_r_down   <- binomial_model(S, K, r - eps_rho, q, sigma, T, N,
                               tipo_opcion, estilo_ejercicio)$precio   # V₋
  
  rho     <- (V_r_up - V_r_down) / (2 * eps_rho)
  rho_pct <- rho / 100   # por 1% de cambio en r
  
  
  return(list(
    delta        = delta,
    gamma        = gamma,
    theta        = theta,
    theta_diario = theta_diario,
    vega         = vega,
    vega_pct     = vega_pct,
    rho          = rho,
    rho_pct      = rho_pct
  ))
}


# GRIEGAS BINOMIALES PARA BARRERAS
# -----------------------------------------------------------------------------
# Idénticas a griegas_binomiales en delta/gamma/theta (diferencias sobre el
# árbol). Para vega y rho se re-ejecuta binomial_barrier_model con perturbación.

griegas_binomiales_barrera <- function(res_bar, S, K, r, q, sigma, T, N,
                                       tipo, estilo, tipo_barrera, B) {
  
  S_tree <- res_bar$S_tree
  V_tree <- res_bar$V_tree
  dt     <- res_bar$dt
  u      <- res_bar$u      # factor de subida
  d      <- res_bar$d      # factor de bajada
  
  V_0     <- V_tree[1, 1]   # V  — valor actual (n=0)
  V_plus  <- V_tree[2, 2]   # V⁺ — valor si precio subió (n=1, j=1)
  V_minus <- V_tree[2, 1]   # V⁻ — valor si precio bajó  (n=1, j=0)
  
  # ---------------------------------------------------------------------------
  # DELTA  
  
  #          V⁺ − V⁻
  #   Δ  =  ─────────
  #          (u − d)·S
  
  # Igual que en vanilla. La diferencia es que V⁺ y V⁻ ya incorporan la
  # condición de barrera pues los nodos muertos tienen V=0 (knock-out), nodos no
  # activados reflejan la probabilidad futura de activación (knock-in).
  # Por eso el delta puede ser negativo cerca de la barrera (por ejemplo up-and-out
  # call -> subir el precio acerca a la barrera y destruye valor).
  # ---------------------------------------------------------------------------
  delta <- (V_plus - V_minus) / ((u - d) * S)
  
  
  # ---------------------------------------------------------------------------
  # GAMMA  
  
  # Γ = ∂²V/∂S² ≈ (Δ⁺ − Δ⁻) / distancia entre los puntos G
  
  # Misma fórmula que vanilla. La gamma puede ser negativa para opciones
  # knock-out cerca de la barrera. La curvatura del valor se invierte porque
  # acercarse a la barrera reduce el valor de forma no lineal.
  # ---------------------------------------------------------------------------
  if (N >= 2) {
    delta_up   <- (V_tree[3, 3] - V_tree[3, 2]) /
      (S_tree[3, 3] - S_tree[3, 2])   # Δ⁺ (nodos superiores de G)
    
    delta_down <- (V_tree[3, 2] - V_tree[3, 1]) /
      (S_tree[3, 2] - S_tree[3, 1])   # Δ⁻ (nodos inferiores de G)
    
    gamma <- (delta_up - delta_down) /
      (0.5 * (S_tree[3, 3] - S_tree[3, 1]))
  } else {
    gamma <- NA
  }
  
  
  # ---------------------------------------------------------------------------
  # THETA
  
  # En CRR, u·d = 1 → el nodo central n=2, j=1 tiene precio exactamente S.
  # Comparamos V al mismo precio pero 2δt más tarde:
  
  #   Θ = (V[n=2, j=1] − V[n=0, j=0]) / (2·δt)
  
  # El nodo central (n=2, j=1) puede ser alcanzado por la ruta subida-bajada
  # o bajada-subida. El árbol barrera ya determina si ese nodo es vivo/muerto
  # (knock-out) o activado/no activado (knock-in), por lo que el theta
  # captura automáticamente el efecto temporal de la barrera.
  # ---------------------------------------------------------------------------
  if (N >= 2) {
    theta        <- (V_tree[3, 2] - V_tree[1, 1]) / (2 * dt)
    theta_diario <- theta / 365   # convertido a €/día
  } else {
    theta        <- NA
    theta_diario <- NA
  }
  
  
  # ---------------------------------------------------------------------------
  # VEGA
  
  #   ν  ≈  (V₊ − V₋) / (2ε)
  
  # Se reconstruye el árbol barrera completo con σ±ε.
  # Nota: para opciones knock-out, la vega puede ser NEGATIVA porque mayor
  # volatilidad aumenta la probabilidad de tocar la barrera y ser eliminada.
  # ---------------------------------------------------------------------------
  eps_vega     <- 0.001    # ε = 0.1% de volatilidad
  V_sigma_up   <- binomial_barrier_model(S, K, r, q, sigma + eps_vega, T, N,
                                         tipo, estilo, tipo_barrera, B)$precio
  V_sigma_down <- binomial_barrier_model(S, K, r, q, sigma - eps_vega, T, N,
                                         tipo, estilo, tipo_barrera, B)$precio
  
  vega     <- (V_sigma_up - V_sigma_down) / (2 * eps_vega)
  vega_pct <- vega / 100   # por 1% de cambio en σ
  
  
  # ---------------------------------------------------------------------------
  # RHO
  
  #   ρ  ≈  (V₊ − V₋) / (2ε)
  
  # Se reconstruye el árbol barrera completo con r±ε.
  # r afecta a la probabilidad riesgo-neutral p y al factor de descuento,
  # y también cambia la geometría del árbol (u, d, p), por eso requiere
  # re-valoración completa.
  # ---------------------------------------------------------------------------
  eps_rho    <- 0.0001    # ε = 0.01%
  V_r_up     <- binomial_barrier_model(S, K, r + eps_rho, q, sigma, T, N,
                                       tipo, estilo, tipo_barrera, B)$precio
  V_r_down   <- binomial_barrier_model(S, K, r - eps_rho, q, sigma, T, N,
                                       tipo, estilo, tipo_barrera, B)$precio
  
  rho     <- (V_r_up - V_r_down) / (2 * eps_rho)
  rho_pct <- rho / 100   # por 1% de cambio en r
  
  
  return(list(
    delta        = delta,
    gamma        = gamma,
    theta        = theta,
    theta_diario = theta_diario,
    vega         = vega,
    vega_pct     = vega_pct,
    rho          = rho,
    rho_pct      = rho_pct
  ))
}


# -----------------------------------------------------------------------------
# 2.3 BLACK-SCHOLES (solución analítica, solo para europeas)
# -----------------------------------------------------------------------------

# (Wilmott, Cap. 7-8):

# La ecuación diferencial de Black-Scholes (1973):
#   ∂V/∂t + ½σ²S²·∂²V/∂S² + rS·∂V/∂S - rV = 0 (Wilmott, ecuación 8.1, p. 170)

# Solución analítica con dividendo continuo q (Merton, 1973; Wilmott, Cap. 8):

#   CALL: C = S·e^(-qT)·N(d₁) - K·e^(-rT)·N(d₂)
#   PUT:  P = K·e^(-rT)·N(-d₂) - S·e^(-qT)·N(-d₁)

# donde:
#   d₁ = [ln(S/K) + (r - q + σ²/2)·T] / (σ·√T)
#   d₂ = d₁ - σ·√T
#   N(·) : función de distribución acumulada de la normal estándar

# INTUICIÓN DE d₁ y d₂ (Wilmott, p. 175-181):
#   d₂ → probabilidad riesgo-neutral de que la opción expire ITM
#   d₁ → d₂ ajustado por la convexidad del payoff (delta de la call)

# Black-Scholes solo es válido para europeas.
# Las americanas requieren métodos numéricos (binomial, diferencias finitas).
# -----------------------------------------------------------------------------

black_scholes <- function(S, K, r, q, sigma, T, tipo_opcion) {
  # Calculamos d₁ y d₂ (wilmott página 177)
  d1 <- (log(S / K) + (r - q + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  if (tipo_opcion == "call") {
    # Fórmula BS para call (Wilmott, p. 177)
    precio <- S * exp(-q * T) * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  } else {  # put
    # Fórmula BS para put (Wilmott, p. 179)
    precio <- K * exp(-r * T) * pnorm(-d2) - S * exp(-q * T) * pnorm(-d1)
  }
  
  return(list(precio = precio, d1 = d1, d2 = d2))
}


# -----------------------------------------------------------------------------
# 2.4 GRIEGAS ANALÍTICAS DE BLACK-SCHOLES
# -----------------------------------------------------------------------------

# (Wilmott, Cap. 8, pp. 182-200):

# Definiciones:
#   N'(x) = φ(x) = (1/√(2π))·e^(-x²/2)   -> densidad normal estándar
#   N(x) : función de distribución acumulada normal estándar

# DELTA (Δ = ∂V/∂S):
#   Call: Δ = e^(-qT)·N(d₁)           [entre 0 y 1]
#   Put:  Δ = e^(-qT)·(N(d₁) - 1)     [entre -1 y 0]
#   → Ratio de cobertura; en el límite q=0: Δ_call = N(d₁)

# GAMMA (Γ = ∂²V/∂S²):
#   Call = Put: Γ = e^(-qT)·N'(d₁) / (S·σ·√T)   (siempre positivo)
#   → Velocidad de cambio del delta; en el límite q=0: Γ = N'(d₁)/(S·σ·√T)

# THETA (Θ = ∂V/∂t):
#   Call: Θ = -S·e^(-qT)·N'(d₁)·σ/(2√T)  +  q·S·e^(-qT)·N(d₁)  -  r·K·e^(-rT)·N(d₂)
#   Put:  Θ = -S·e^(-qT)·N'(d₁)·σ/(2√T)  -  q·S·e^(-qT)·N(-d₁) +  r·K·e^(-rT)·N(-d₂)
#   → Decaimiento temporal (normalmente negativo para posiciones largas)

# VEGA (ν = ∂V/∂σ):
#   Call = Put: ν = S·e^(-qT)·√T·N'(d₁)   (siempre positivo)
#   → Sensibilidad a la volatilidad

# RHO (ρ = ∂V/∂r):
#   Call: ρ = K·T·e^(-rT)·N(d₂)
#   Put:  ρ = -K·T·e^(-rT)·N(-d₂)
#   → Sensibilidad al tipo de interés libre de riesgo
# -----------------------------------------------------------------------------

griegas_bs <- function(S, K, r, q, sigma, T, tipo_opcion) {
  d1 <- (log(S / K) + (r - q + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  # Densidad normal estándar evaluada en d₁: N'(d₁) = φ(d₁)
  phi_d1 <- dnorm(d1)   # equivalente a (1/√(2π))·exp(-d1²/2)
  
  # DELTA
  if (tipo_opcion == "call") {
    delta <- exp(-q * T) * pnorm(d1)           # e^(-qT) N(d₁)
  } else {
    delta <- exp(-q * T) * (pnorm(d1) - 1)     # e^(-qT) (N(d₁) - 1)
  }
  
  # GAMMA (igual para call y put)
  gamma <- exp(-q * T) * phi_d1 / (S * sigma * sqrt(T))
  
  # THETA 
  if (tipo_opcion == "call") {
    theta <- (-S * exp(-q * T) * phi_d1 * sigma / (2 * sqrt(T))) +
      q * S * exp(-q * T) * pnorm(d1) -
      r * K * exp(-r * T) * pnorm(d2)
  } else {
    theta <- (-S * exp(-q * T) * phi_d1 * sigma / (2 * sqrt(T))) -
      q * S * exp(-q * T) * pnorm(-d1) +
      r * K * exp(-r * T) * pnorm(-d2)
  }
  theta_diario <- theta / 365
  
  # VEGA (igual para call y put, expresada por unidad de σ)
  vega <- S * exp(-q * T) * sqrt(T) * phi_d1
  vega_pct <- vega / 100  # por 1% de cambio en σ
  
  # RHO
  if (tipo_opcion == "call") {
    rho <- K * T * exp(-r * T) * pnorm(d2)
  } else {
    rho <- -K * T * exp(-r * T) * pnorm(-d2)
  }
  rho_pct <- rho / 100
  
  return(list(
    delta        = delta,
    gamma        = gamma,
    theta        = theta,
    theta_diario = theta_diario,
    vega         = vega,
    vega_pct     = vega_pct,
    rho          = rho,
    rho_pct      = rho_pct,
    d1           = d1,
    d2           = d2
  ))
}


# -----------------------------------------------------------------------------
# MÉTODO DE AJUSTE DEL SPOT INICIAL PARA DIVIDENDOS DISCRETOS
# -----------------------------------------------------------------------------

# Haug, Haug & Lewis, 2003; Hull, 2000:
# Los dividendos discretos son pagos puntuales D_i en fechas t_i que reducen
# el precio de la acción en la misma cantidad cuando se pagan. 
# Una aproximación válida para opciones europeas vanilla consiste en restar al precio 
# inicial S₀ el valor presente de todos los dividendos futuros dentro de la vida de la opción:

#   S* = S₀ - Σᵢ Dᵢ · exp(-r · tᵢ)

# El árbol binomial (o Black-Scholes / Monte Carlo) se construye con S* como
# precio inicial, y la valoración procede de forma estándar.
# Al no alterar la geometría del árbol, sigue siendo recombinante.

# Hay que tener en cuenta que:
#   - Solo es fiable para opciones europeas vanilla (no americanas, no
#     path-dependent), ya que altera la estructura temporal de precios.
#   - Funciona mejor con dividendos pequeños y vencimientos cortos.


# REFERENCIA: Haug, Haug & Lewis (2003) "Back to Basics: a new approach to
# the discrete dividend problem" — el método M73 (escrowed dividend) es
# equivalente a este ajuste de spot bajo el modelo GBM.
# -----------------------------------------------------------------------------

# Calcula el spot ajustado S* restando el valor presente de los dividendos discretos.
# Entradas:
#   S      : precio spot actual S₀
#   r      : tasa libre de riesgo anualizada (continua)
#   divs   : data.frame con columnas 'importe' (D_i) y 'tiempo' (t_i en años)
#            Solo se consideran dividendos con 0 < t_i <= T
#   T      : tiempo hasta vencimiento (años)
# Salida: S* (precio spot ajustado)

calcular_spot_ajustado <- function(S, r, divs, T) {
  
  divs_validos <- divs[divs$tiempo > 0 & divs$tiempo <= T, ]
  
  if (nrow(divs_validos) == 0) return(S)  # sin dividendos: S* = S₀
  
  # Suma de valores presentes: Σ Dᵢ · exp(-r · tᵢ)
  pv_divs <- sum(divs_validos$importe * exp(-r * divs_validos$tiempo))
  
  # S* = S₀ - PV(dividendos); se garantiza S* > 0
  S_star <- S - pv_divs
  return(S_star)  
}

# Black-Scholes con ajuste de spot para dividendos discretos (europeas vanilla).
# Se calcula S* = S₀ - Σ Dᵢ exp(-r tᵢ) y se aplica la fórmula BS estándar con S*.

black_scholes_spot_adj <- function(S, K, r, sigma, T, tipo_opcion, divs) {
  # Con ajuste de spot no se usa tasa de dividendos q (q=0 en BS):
  # el efecto de los dividendos ya está incorporado en S*.
  S_star <- calcular_spot_ajustado(S, r, divs, T)
  
  # Aplicamos Black-Scholes estándar con q=0 y spot ajustado S*
  # (Haug, Haug & Lewis, 2003, ecuación M73)
  resultado <- black_scholes(S_star, K, r, q = 0, sigma, T, tipo_opcion)
  resultado$S_star <- S_star   # devolvemos S* para mostrarlo en la interfaz
  return(resultado)
}

# Árbol binomial con ajuste de spot para dividendos discretos (europeas vanilla).
# Se construye el árbol con S* en lugar de S₀; el resto es idéntico al binomial estándar.

binomial_model_spot_adj <- function(S, K, r, sigma, T, N, tipo_opcion, estilo_ejercicio, divs) {
  S_star <- calcular_spot_ajustado(S, r, divs, T)
  # Árbol estándar con S* y q=0 (los dividendos ya están en S*)
  resultado <- binomial_model(S_star, K, r, q = 0, sigma, T, N, tipo_opcion, estilo_ejercicio)
  resultado$S_star <- S_star
  return(resultado)
}

# Monte Carlo con ajuste de spot para dividendos discretos (europeas vanilla).
# Se simula con S* y q=0; el payoff se calcula igual que en mc_vanilla.

mc_vanilla_spot_adj <- function(S, K, r, sigma, T, n_sims, tipo_opcion, semilla, divs) {
  S_star <- calcular_spot_ajustado(S, r, divs, T)
  # mc_vanilla estándar con spot ajustado y q=0
  resultado <- mc_vanilla(S_star, K, r, q = 0, sigma, T, n_sims, tipo_opcion, semilla)
  resultado$S_star <- S_star
  return(resultado)
}


# -----------------------------------------------------------------------------
# 2.5 FUNCIÓN PARA PREPARAR DATOS DEL ÁRBOL PARA GGPLOT
# -----------------------------------------------------------------------------
# Convierte las matrices S_tree y V_tree en un data.frame largo para ggplot2.
# Cada fila representa un nodo del árbol con su posición (paso, nivel)
# y sus valores de S y V.

preparar_datos_arbol <- function(resultado, K, tipo_opcion = "call", estilo_ejercicio = "european") {
  S_tree <- resultado$S_tree
  V_tree <- resultado$V_tree
  N      <- resultado$N
  intrinsic_tree <- matrix(0, nrow = N + 1, ncol = N + 1)
  continuation_tree <- matrix(NA, nrow = N + 1, ncol = N + 1)
  exercise_tree <- matrix(FALSE, nrow = N + 1, ncol = N + 1)
  
  if (!is.null(resultado$intrinsic_tree)) {
    intrinsic_tree <- resultado$intrinsic_tree
  } else {
    for (n in 0:N) {
      for (j in 0:n) {
        S_val <- S_tree[n + 1, j + 1]
        intrinsic_tree[n + 1, j + 1] <- if (tipo_opcion == "call") max(S_val - K, 0) else max(K - S_val, 0)
      }
    }
  }
  
  if (!is.null(resultado$continuation_tree)) {
    continuation_tree <- resultado$continuation_tree
  }
  
  if (!is.null(resultado$exercise_tree)) {
    exercise_tree <- resultado$exercise_tree
  }
  
  show_exercise_info <- estilo_ejercicio == "american"
  
  nodos <- data.frame()
  for (n in 0:N) {
    for (j in 0:n) {
      S_val <- S_tree[n + 1, j + 1]
      V_val <- V_tree[n + 1, j + 1]
      intrinsic <- intrinsic_tree[n + 1, j + 1]
      continuation <- continuation_tree[n + 1, j + 1]
      itm <- intrinsic > 0
      early <- show_exercise_info && isTRUE(exercise_tree[n + 1, j + 1])
      
      hover_base <- paste0(
        "<br>Underlying price: ", round(S_val, 4),
        "<br>Option value: ", round(V_val, 6),
        "<br>Status: ", ifelse(early, "Early exercise", ifelse(itm, "ITM", "OTM"))
      )
      
      hover_exercise <- if (show_exercise_info) {
        paste0(
          "<br>Intrinsic value: ", round(intrinsic, 6),
          "<br>Continuation value: ", ifelse(is.na(continuation), "N/A", round(continuation, 6)),
          "<br>Early exercise: ", ifelse(early, "Yes", "No")
        )
      } else {
        ""
      }
      
      nodos <- rbind(nodos, data.frame(
        paso  = n,
        nivel = j,
        S     = round(S_val, 4),
        V     = round(V_val, 6),
        intrinsic = round(intrinsic, 6),
        continuation = ifelse(is.na(continuation), NA, round(continuation, 6)),
        ITM   = itm,
        early_exercise = early,
        y_pos = j - n / 2,
        hover = paste0(hover_base, hover_exercise),
        stringsAsFactors = FALSE
      ))
    }
  }
  nodos
}

# PREPARAR DATOS DEL ÁRBOL BARRERA PARA GGPLOT
# -----------------------------------------------------------------------------

preparar_datos_arbol_barrera <- function(res_bar, K, tipo_opcion = "call", estilo_ejercicio = "european") {
  S_tree <- res_bar$S_tree
  V_tree <- res_bar$V_tree
  N      <- res_bar$N
  alive  <- res_bar$alive_matrix
  es_ko  <- res_bar$es_ko
  exercise_tree <- res_bar$exercise_tree
  intrinsic_tree <- res_bar$intrinsic_tree
  continuation_tree <- if (!is.null(res_bar$continuation_tree)) res_bar$continuation_tree else matrix(NA, nrow = N + 1, ncol = N + 1)
  show_exercise_info <- estilo_ejercicio == "american"
  
  nodos <- data.frame()
  for (n in 0:N) {
    for (j in 0:n) {
      active_state <- alive[n + 1, j + 1]
      estado <- if (es_ko) {
        if (active_state) "Alive" else "Knocked out"
      } else {
        if (active_state) "Activated" else "Not activated"
      }
      S_val <- S_tree[n + 1, j + 1]
      V_val <- V_tree[n + 1, j + 1]
      intrinsic <- intrinsic_tree[n + 1, j + 1]
      continuation <- continuation_tree[n + 1, j + 1]
      # El ejercicio anticipado solo tiene sentido en nodos donde la opción está realmente activa.
      # Para knock-out: el nodo debe seguir vivo; los nodos eliminados valen cero y no pueden ejercerse.
      # Para knock-in:  el nodo debe haber sido ya activado; los nodos sin activar no pueden ejercerse aún.
      early <- show_exercise_info && isTRUE(exercise_tree[n + 1, j + 1]) && isTRUE(active_state)
      v_active <- if (!is.null(res_bar$V_active)) res_bar$V_active[n + 1, j + 1] else NA
      v_not <- if (!is.null(res_bar$V_not_active)) res_bar$V_not_active[n + 1, j + 1] else NA
      
      extra_state <- if (es_ko) {
        paste0("<br>Path state: ", estado)
      } else {
        paste0("<br>Path state: ", estado,
               "<br>Value if activated: ", ifelse(is.na(v_active), "N/A", round(v_active, 6)),
               "<br>Value if not activated: ", ifelse(is.na(v_not), "N/A", round(v_not, 6)))
      }
      
      hover_base <- paste0(
        "<br>Underlying price: ", round(S_val, 4),
        "<br>Option value: ", round(V_val, 6),
        extra_state
      )
      
      hover_exercise <- if (show_exercise_info) {
        paste0(
          "<br>Intrinsic value: ", round(intrinsic, 6),
          "<br>Continuation value: ", ifelse(is.na(continuation), "N/A", round(continuation, 6)),
          "<br>Early exercise: ", ifelse(early, "Yes", "No")
        )
      } else {
        ""
      }
      
      nodos <- rbind(nodos, data.frame(
        paso = n,
        nivel = j,
        S = round(S_val, 4),
        V = round(V_val, 6),
        intrinsic = round(intrinsic, 6),
        continuation = ifelse(is.na(continuation), NA, round(continuation, 6)),
        active = active_state,
        estado = estado,
        early_exercise = early,
        y_pos = j - n / 2,
        hover = paste0(hover_base, hover_exercise),
        stringsAsFactors = FALSE
      ))
    }
  }
  nodos
}

# -----------------------------------------------------------------------------
# 2.6 FUNCIONES MONTE CARLO
# -----------------------------------------------------------------------------

# FÓRMULA FUNDAMENTAL (Wilmott, p. 148):

#   V(S, t) = e^{-r(T-t)} · E*[ Payoff(S(T)) ]

#   El precio de cualquier derivado es el valor presente del payoff eesperado
#   bajo la medida riesgo-neutral (probabilidades p*).
#   Monte Carlo aproxima esa esperanza con la media de M simulaciones:

#     V ≈ e^{-rT} · (1/M) · Σᵢ Payoff(Sᵢ(T))

#   La clave es que los caminos se generan con la dinámica riesgo-neutral:
#     dS = (r-q)·S·dt + σ·S·dX    (drift riesgo-neutral con dividendos)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 2.6.1 MÉTODO EXACTO (Wilmott, ec. 29.1, p. 585)
# -----------------------------------------------------------------------------

# La ecuación diferencial estocástica riesgo-neutral:
#   dS = (r-q)·S·dt + σ·S·dX

# Reescrita en términos de log(S):
#   d(log S) = (r - q - σ²/2)·dt + σ·dX

# Integrando exactamente sobre un paso δt (sin aproximación):

#   S(t + δt) = S(t) · exp[ (r - q - σ²/2)·δt  +  σ·√δt·φ ]     (ec. 29.1)

#   donde φ ~ N(0,1) es una variable aleatoria normal estándar.

# Esta expresión es exacta para cualquier δt.
# Para opciones vanilla (solo dependen de S(T)), podemos usar δt = T
# y dar solo salto desde t=0 hasta T, sin error de discretización.
# Para opciones path-dependent necesitamos pasos pequeños (por ejemplo δt = T/252).

# Comparamos con el método de Euler (aproximación):
#   S(t + δt) ≈ S(t) + (r-q)·S(t)·δt + σ·S(t)·√δt·φ    
#   → Error O(δt); menos preciso para el mismo δt.
#   Usamos los mismos números φ al comparar ambos métodos.
## Sólo comparo con el método Euler para opciones vanilla europeas, tiene sentido hacerlo para algún otro caso?
# -----------------------------------------------------------------------------

mc_simular_gbm <- function(S, r, q, sigma, T, n_steps, n_sims, semilla = NULL) {
  # Genera una matriz [n_sims × (n_steps+1)] de trayectorias de precios
  # usando la ecuación 29.1 de Wilmott.
  
  # S: precio inicial S₀
  # r: tasa libre de riesgo
  # q: rentabilidad por dividendo continua anualizada
  # sigma: volatilidad σ
  # T: tiempo hasta vencimiento (años)
  # n_steps: número de pasos temporales
  #     → Para vanilla: usa n_steps = 1 (un salto exacto de δt = T)
  #     → Para path-dependent: usa n_steps = número indicado por el usuario 
  # n_sims: número de simulaciones (caminos)
  # semilla: para reproducibilidad 
  
  # Comprobamos si se ha proporcionado una semilla aleatoria
  # !is.null(semilla) verifica que el valor de semilla existe
  # semilla > 0 indica que el usuario quiere fijar una semilla concreta
  # Si ambas condiciones se cumplen, set.seed fija la semilla
  # Esto permite que la simulación Monte Carlo sea reproducible
  if (!is.null(semilla) && semilla > 0) set.seed(semilla)
  
  dt <- T / n_steps
  
  # Coeficientes de la ec. 29.1:
  drift <- (r - q - 0.5 * sigma^2) * dt   # término determinista
  difusion <- sigma * sqrt(dt)            # coeficiente del término estocástico
  
  # Matriz de números aleatorios normales estándar
  # Dimensión: n_sims filas × n_steps columnas
  phi <- matrix(rnorm(n_sims * n_steps), nrow = n_sims, ncol = n_steps)
  
  # Construir trayectorias aplicando ec. 29.1 en cada paso
  caminos <- matrix(0, nrow = n_sims, ncol = n_steps + 1)
  caminos[, 1] <- S   # todos los caminos empiezan en S₀
  
  for (paso in 1:n_steps) {
    # S(t + δt) = S(t) · exp[ (r - q - σ²/2)·δt + σ·√δt·φ ]
    caminos[, paso + 1] <- caminos[, paso] *
      exp(drift + difusion * phi[, paso])
  }
  
  # Devolvemos caminos y la matriz phi para poder reusar los mismos números
  # en la comparación con el método de Euler
  return(list(caminos = caminos, phi = phi, dt = dt))
}


# MC OPCIONES VANILLA (Call / Put)

# Las opciones vanilla son path-independent pues el payoff solo depende de S(T).
# Por tanto podemos dar un solo salto exacto con δt = T (Wilmott, p. 585):
#   "if we have a payoff that only depends on the final asset value...
#    we can simulate the final asset price in one giant leap"

# Payoffs:
#   Call: max(S(T) - K, 0)
#   Put:  max(K - S(T), 0)

# Precio MC: V ≈ e^{-rT} · mean(payoffs)

mc_vanilla <- function(S, K, r, q, sigma, T, n_sims, tipo, semilla = NULL) {
  # Un solo paso (δt = T): exacto, sin error de discretización
  res <- mc_simular_gbm(S, r, q, sigma, T, n_steps = 1, n_sims, semilla)
  ST  <- res$caminos[, 2]   # precios en T (columna 2 de la matriz 2-columnas)
  
  payoffs <- if (tipo == "call") pmax(ST - K, 0) else pmax(K - ST, 0)
  
  # Devolvemos payoffs para calcular IC y caminos para visualización
  list(payoffs = payoffs, ST = ST, caminos = res$caminos, phi = res$phi)
}


# MC OPCIONES ASIÁTICAS

# Strongly path-dependent (Wilmott, p. 265): el payoff depende del precio
# medio del subyacente durante la vida de la opción, no solo del precio final.
# Necesitamos simular la trayectoria completa.

# DOS TIPOS DE MEDIA:

#   Media Aritmética:
#      A = (1/n) · Σ S(tᵢ)
#      Payoff Call: max(A - K, 0)
#      Payoff Put:  max(K - A, 0)
#      → No existe fórmula cerrada exacta.
#        Aproximación analítica: método M1/M2

#   Media Geométrica (Wilmott, p. 268):
#      G = exp[ (1/n) · Σ log S(tᵢ) ] 
#      → La distribución de G también es lognormal (volatilidad reducida σ/√3)
#      → Sí existe fórmula cerrada exacta 

mc_asiatica <- function(S, K, r, q, sigma, T, n_steps, n_sims, tipo,
                        tipo_media, semilla = NULL) {
  
  res     <- mc_simular_gbm(S, r, q, sigma, T, n_steps, n_sims, semilla)
  caminos <- res$caminos   # [n_sims × (n_steps+1)]
  
  # Calcular la media de cada camino (fila)
  if (tipo_media == "aritmetica") {
    # Media aritmética de todos los precios del camino (t=0 incluido)
    medias <- rowMeans(caminos)
  } else {
    # Media geométrica: exp(media de log-precios)  - Wilmott, p. 268
    # La media geométrica es exacta para distribuciones lognormales
    medias <- exp(rowMeans(log(caminos)))
  }
  
  payoffs <- if (tipo == "call") pmax(medias - K, 0) else pmax(K - medias, 0)
  
  list(payoffs = payoffs, caminos = caminos, medias = medias)
}


# MC OPCIONES BARRERA

# Weakly path-dependent (Wilmott, p. 264): el payoff depende de si el precio
# ha cruzado un nivel barrera B en algún momento de la vida de la opción.

# TIPOS:
#   Knock-Out: opción desaparece si toca la barrera → payoff = 0
#     Down-and-Out: B < S₀ (barrera por debajo)
#     Up-and-Out:   B > S₀ (barrera por encima)
#   Knock-In: opción nace si toca la barrera → payoff = vanilla si toca, 0 si no
#     Down-and-In:  B < S₀
#     Up-and-In:    B > S₀

# PARIDAD: Knock-In + Knock-Out = Vanilla (misma B, mismo K)
# Se puede usar como comprobación de los resultados.

# LIMITACIÓN DE MC: con pasos discretos podemos no detectar que el precio 
# cruzó la barrera entre dos pasos. Error O(√δt).
# La solución es usar muchos pasos (n_steps = 252 o más).

mc_barrera <- function(S, K, r, q, sigma, T, n_steps, n_sims, tipo,
                       tipo_barrera, nivel_barrera, semilla = NULL) {
  # tipo_barrera: "down-and-out", "down-and-in", "up-and-out", "up-and-in"
  # nivel_barrera: B (nivel de la barrera)
  
  res     <- mc_simular_gbm(S, r, q, sigma, T, n_steps, n_sims, semilla)
  caminos <- res$caminos
  B       <- nivel_barrera
  
  # Comprobar si cada camino tocó la barrera (excluimos t=0)
  caminos_post_inicio <- caminos[, -1]
  
  toco_down <- apply(caminos_post_inicio, 1, function(path) any(path <= B)) 
  toco_up   <- apply(caminos_post_inicio, 1, function(path) any(path >= B))
  
  # Precio final en T
  ST <- caminos[, n_steps + 1]
  
  # Payoff vanilla base (sin condición de barrera)
  payoff_base <- if (tipo == "call") pmax(ST - K, 0) else pmax(K - ST, 0)
  
  # Aplicar condición de barrera
  payoffs <- switch(tipo_barrera,
                    "down-and-out" = ifelse(toco_down, 0,           payoff_base),
                    "down-and-in"  = ifelse(toco_down, payoff_base, 0),
                    "up-and-out"   = ifelse(toco_up,   0,           payoff_base),
                    "up-and-in"    = ifelse(toco_up,   payoff_base, 0)
  )
  
  list(payoffs = payoffs, caminos = caminos,
       toco = if (grepl("down", tipo_barrera)) toco_down else toco_up)
}


# MC OPCIONES LOOKBACK

# Strongly path-dependent (Wilmott, p. 268): el payoff depende del máximo
# o mínimo realizado del precio durante toda la vida de la opción. 

# DOS VARIANTES:

#   Floating Strike (el strike es el extremo del camino):
#      Call: S(T) - min_path(S)     [siempre ≥ 0, el min siempre ≤ S(T)]
#      Put:  max_path(S) - S(T)     [siempre ≥ 0, el max siempre ≥ S(T)]
#      → Existe fórmula cerrada (Wilmott, pp. 268-269) 

#   Fixed Strike K (el extremo reemplaza a S(T) en el payoff):
#      Call: max(max_path(S) - K, 0)   → se compara el máximo alcanzado con K
#      Put:  max(K - min_path(S), 0)   → se compara K con el mínimo alcanzado
#      → Existe fórmula cerrada (Wilmott, p. 269) 

mc_lookback <- function(S, K, r, q, sigma, T, n_steps, n_sims, tipo,
                        tipo_strike, semilla = NULL) {
  # tipo_strike: "floating" o "fixed"
  
  res     <- mc_simular_gbm(S, r, q, sigma, T, n_steps, n_sims, semilla)
  caminos <- res$caminos
  
  ST         <- caminos[, n_steps + 1]
  min_camino <- apply(caminos, 1, min)
  max_camino <- apply(caminos, 1, max)
  
  if (tipo_strike == "floating") {
    
    if (tipo == "call") {
      payoffs <- pmax(ST - min_camino, 0)
    } else {
      payoffs <- pmax(max_camino - ST, 0)}
  } else {
    if (tipo == "call") {
      payoffs <- pmax(max_camino - K, 0)
    } else {
      payoffs <- pmax(K - min_camino, 0)}}
  
  list(payoffs = payoffs, caminos = caminos,
       min_camino = min_camino, max_camino = max_camino)
}


# -----------------------------------------------------------------------------
# 2.6.2 FÓRMULAS CERRADAS (benchmark para comparar con MC)
# -----------------------------------------------------------------------------

# --- Asiática Geométrica (Wilmott, p. 267) ---

# "The geometric average of a lognormal random walk is itself lognormally
#  distributed, but with a reduced volatility."

# Al inicio (t=0, sin historia previa):
#   I = ∫₀⁰ log S(τ)dτ = 0  →  G = S (el "promedio hasta ahora" es S₀)

# Denominador: σ·√(T³/3)   [Wilmott: σ√((T-t)³/3) evaluado en t=0]
# d₁ = [ T·log(S/K) + (r-q-σ²/2)·T²/2 + σ²·T²/3 ] / (σ·√(T³/3))
# d₂ = [ T·log(S/K) + (r-q-σ²/2)·T²/2           ] / (σ·√(T³/3))

# Call = e^{-rT}·[ S·exp((r -q -σ²/2)·T/2 + σ²·T/6)·N(d₁) - K·N(d₂) ]
# Put  = e^{-rT}·[ K·N(-d₂) - S·exp(...)·N(-d₁) ]

bs_asiatica_geometrica <- function(S, K, r, q, sigma, T) {
  # Factor exponencial del "promedio" esperado al inicio
  exp_adj <- exp((r - q - sigma^2 / 2) * T / 2 + sigma^2 * T / 6)
  
  denom <- sigma * sqrt(T^3 / 3)   # denominador de d₁ y d₂
  
  d1 <- (T * log(S / K) + (r - q - sigma^2 / 2) * T^2 / 2 + sigma^2 * T^2 / 3) / denom
  d2 <- (T * log(S / K) + (r - q - sigma^2 / 2) * T^2 / 2                     ) / denom
  
  call <- exp(-r * T) * (S * exp_adj * pnorm(d1) - K * pnorm(d2))
  put  <- exp(-r * T) * (K * pnorm(-d2) - S * exp_adj * pnorm(-d1)) ## en la página 268 de Wilmott pone d1 en lugar de -d1, entiendo que es un error del libro?
  
  list(call = call, put = put, d1 = d1, d2 = d2)
}

# --- Asiática Aritmética: aproximación M1/M2 --- 

## https://www-2.rotman.utoronto.ca/~hull/technicalnotes/TechnicalNote27.pdf
## https://homepage.ntu.edu.tw/~jryanwang/courses/Financial%20Computation%20or%20Financial%20Engineering%20%28graduate%20level%29/FE_Ch10%20Asian%20Option.pdf
## https://dspacemainprd01.lib.uwaterloo.ca/server/api/core/bitstreams/1f7397de-e49f-4d16-ad8a-89d080dc352c/content

# Para las opciones asiáticas de media aritmética no existe, en general,
# una fórmula cerrada exacta bajo el modelo de Black-Scholes.

# A diferencia de la media geométrica, la media aritmética de precios lognormales
# no sigue una distribución lognormal exacta. Por ello, una alternativa habitual
# consiste en utilizar una aproximación por ajuste de momentos.

# La idea de esta aproximación es:
#   1. Calcular el primer momento del promedio aritmético, M1 = E[A_T].
#   2. Calcular el segundo momento del promedio aritmético, M2 = E[A_T^2].
#   3. Aproximar la distribución del promedio aritmético mediante una lognormal
#      con los mismos dos primeros momentos.
#   4. Valorar la opción mediante una expresión análoga a Black-Scholes,
#      utilizando una volatilidad ajustada.

# En este caso se trabaja con un promedio aritmético continuo:

#  A_T = (1/T) ∫_0^T S_t dt

# Además, se permite la existencia de dividendos continuos mediante:

#  b = r - q

# donde:
#   S = precio inicial del activo subyacente
#   K = precio de ejercicio
#   r = tipo de interés libre de riesgo
#   q = rentabilidad por dividendo continua
#   sigma = volatilidad del subyacente
#   T = tiempo hasta vencimiento, expresado en años

bs_asiatica_aritmetica <- function(S, K, r, q, sigma, T) {
  
  b <- r - q
  
  # Función auxiliar:
  # calcula (exp(aT) - 1) / a de forma estable.
  
  # Se utiliza expm1(a * T) en lugar de exp(a * T) - 1 porque es más preciso
  # cuando a está cerca de cero.
  
  # Si a es prácticamente cero, se aplica el límite:
  #  lim_{a -> 0} (e^{aT} - 1) / a = T
  
  expm1_over_a <- function(a) {
    if (abs(a) < 1e-12) {
      return(T)
    } else {
      return(expm1(a * T) / a)
    }
  }
  
  # Primer momento del promedio aritmético continuo.
  #   M1 = E[A_T]
  
  # Para b distinto de cero:
  #  M1 = S * (e^{bT} - 1) / (bT)
  
  # Si b es aproximadamente cero, se aplica el límite y M1 = S.
  
  if (abs(b) < 1e-12) {
    M1 <- S
  } else {
    M1 <- S * (exp(b * T) - 1) / (b * T)
  }
  
  
  # Segundo momento del promedio aritmético continuo.
  #        M2 = E[A_T^2]
  
  # Esta expresión permite obtener la varianza implícita del promedio
  # aritmético y, posteriormente, ajustar una distribución lognormal
  # equivalente.
  
  M2 <- (2 * S^2 / (T^2 * (b + sigma^2))) *
    (expm1_over_a(2 * b + sigma^2) - expm1_over_a(b))
  
  
  # Volatilidad ajustada de la lognormal equivalente.
  
  # Si A_T se aproxima mediante una variable lognormal, la volatilidad
  # equivalente se obtiene igualando los dos primeros momentos:
  
  # sigma_adj = sqrt((1/T) * log(M2 / M1^2))
  
  sig_adj <- sqrt((1 / T) * log(M2 / M1^2))
  
  
  # Si la volatilidad ajustada es prácticamente cero, la opción se valora
  # como el valor actual del payoff determinista sobre el promedio esperado.
  if (sig_adj < 1e-12) {
    
    call <- exp(-r * T) * max(M1 - K, 0)
    put  <- exp(-r * T) * max(K - M1, 0)
    
  } else {
    
    # Parámetros d1 y d2 de la fórmula tipo Black-Scholes.
    
    # M1 se trata como la esperanza del promedio aritmético en vencimiento. 
    # Por eso el payoff completo se descuenta con exp(-rT).
    
    d1 <- (log(M1 / K) + 0.5 * sig_adj^2 * T) /
      (sig_adj * sqrt(T))
    
    d2 <- d1 - sig_adj * sqrt(T)
    
    # Fórmula aproximada para la call asiática aritmética:
    #  C ≈ e^{-rT} [ M1 N(d1) - K N(d2) ]
    
    call <- exp(-r * T) * (M1 * pnorm(d1) - K * pnorm(d2))
    
    
    # Fórmula aproximada para la put asiática aritmética:
    #  P ≈ e^{-rT} [ K N(-d2) - M1 N(-d1) ]
    
    put <- exp(-r * T) * (K * pnorm(-d2) - M1 * pnorm(-d1))
  }
  
  list(
    call = call,
    put = put,
    M1 = M1,
    M2 = M2,
    sigma_adj = sig_adj
  )
}


# --- Lookback Floating Strike (Wilmott, p. 268-269) ---

# Evaluación en t = 0 (inicio de la opción, sin historia previa):
#   M_min = M_max = S₀   →   log(S/M) = 0  cuando d₁ usa M

# Restricción: b = r − q ≠ 0. Si b ≈ 0 (r ≈ q) la corrección
#   σ²/(2b) diverge y la función devuelve NA.

# El strike es el propio extremo realizado del camino:

#   CALL: payoff = max(S(T) − M_min, 0) = S(T) − M_min  [compra al precio mínimo]
#   PUT:  payoff = max(M_max − S(T), 0) = M_max − S(T)  [vende al precio máximo]

# Ambas opciones tienen payoff no negativo en vencimiento.

# Fórmulas de Wilmott en t = 0, D = q, M = S (extremo = precio inicial):

#   d₁ = [log(S/M) + (b + σ²/2)·T] / (σ·√T)
#      = (b + σ²/2)·√T / σ           [log(S/S) = 0]

#   d₂ = d₁ − σ·√T

# CALL (Wilmott, p. 268):
#   V = S·e^{-qT}·N(d₁) − S·e^{-rT}·N(d₂)
#       + S·e^{-rT}·σ²/(2b)·[ N(−d₁ + 2b√T/σ) − e^{bT}·N(−d₁) ]

#   [(S/M)^{-2b/σ²} = 1 ya que M = S]

# PUT (Wilmott, p. 269):
#   V = S·e^{-rT}·N(−d₂) − S·e^{-qT}·N(−d₁)
#       + S·e^{-rT}·σ²/(2b)·[ −N(d₁ − 2b√T/σ) + e^{bT}·N(d₁) ]

bs_lookback_floating <- function(S, r, q, sigma, T) {
  
  b <- r - q   
  
  # Requiere b ≠ 0
  if (abs(b) < 1e-10) {
    message("r ≈ q → b ≈ 0, formula not applicable.")
    return(list(call = NA, put = NA))
  }
  
  # d₁ con M = S → log(S/M) = 0
  d1 <- (b + sigma^2 / 2) * sqrt(T) / sigma
  d2 <- d1 - sigma * sqrt(T)
  
  # Factor común del término de corrección: S·e^{-rT}·σ²/(2b)
  cf <- S * exp(-r * T) * sigma^2 / (2 * b)
  
  call <- S * exp(-q * T) * pnorm(d1) - S * exp(-r * T) * pnorm(d2) +
    cf * (pnorm(-d1 + 2 * b * sqrt(T) / sigma) - exp(b * T) * pnorm(-d1))
  
  put <- S * exp(-r * T) * pnorm(-d2) - S * exp(-q * T) * pnorm(-d1) +
    cf * (-pnorm(d1 - 2 * b * sqrt(T) / sigma) + exp(b * T) * pnorm(d1))
  
  list(call = max(call, 0), put = max(put, 0))
}


# --- Lookback Fixed Strike (Wilmott, p. 269-270) ---

# El strike K es fijo; el payoff depende del extremo histórico:

#   CALL: payoff = max(M_max − K, 0)  [máximo histórico vs strike fijo]
#   PUT:  payoff = max(K − M_min, 0)  [strike fijo vs mínimo histórico]

# Al inicio t = 0: M_max = M_min = S (único precio observado).
# Esto genera dos casos para cada opción según K ≷ S.

# FIXED CALL (M_max = S al inicio): 

#   Caso K ≥ S  [E ≥ M, Wilmott Image 2 "For E > M"]
#     d₁ basado en log(S/K):
#     d₁ = [log(S/K) + (b + σ²/2)·T] / (σ·√T)
#     V = S·e^{-qT}·N(d₁) − K·e^{-rT}·N(d₂)
#         + S·e^{-rT}·σ²/(2b)·[ −(S/K)^{-2b/σ²}·N(d₁ − 2b√T/σ) + e^{bT}·N(d₁) ]

#   Caso K < S  [E < M, Wilmott Image 2 "For E < M"]
#     d₁ basado en log(S/M) = 0:
#     d₁ = (b + σ²/2)·√T / σ
#     V = (S − K)·e^{-rT}                          
#         + S·e^{-qT}·N(d₁) − S·e^{-rT}·N(d₂)
#         + S·e^{-rT}·σ²/(2b)·[ −N(d₁ − 2b√T/σ) + e^{bT}·N(d₁) ]

# FIXED PUT (M_min = S al inicio): 

#   Caso K ≤ S  [E ≤ M, Wilmott Image 2 "For E < M"]
#     d₁ basado en log(S/K):
#     d₁ = [log(S/K) + (b + σ²/2)·T] / (σ·√T)
#     V = K·e^{-rT}·N(−d₂) − S·e^{-qT}·N(−d₁)
#         + S·e^{-rT}·σ²/(2b)·[ (S/K)^{-2b/σ²}·N(−d₁ + 2b√T/σ) − e^{bT}·N(−d₁) ]

#   Caso K > S  [E > M, Wilmott p. 270 "When E > M"]
#     d₁ basado en log(S/M) = 0:
#     d₁ = (b + σ²/2)·√T / σ
#     V = (K − S)·e^{-rT}                         
#         + S·e^{-rT}·N(−d₂) − S·e^{-qT}·N(−d₁)
#         + S·e^{-rT}·σ²/(2b)·[ N(−d₁ + 2b√T/σ) − e^{bT}·N(−d₁) ]

bs_lookback_fixed <- function(S, K, r, q, sigma, T) {
  
  b <- r - q
  
  if (abs(b) < 1e-10) {
    message("r ≈ q → b ≈ 0, formula not applicable.")
    return(list(call = NA, put = NA))
  }
  
  # d₁ y d₂ basados en log(S/K) ]
  d1K <- (log(S / K) + (b + sigma^2 / 2) * T) / (sigma * sqrt(T))
  d2K <- d1K - sigma * sqrt(T)
  
  # d₁ y d₂ basados en log(S/S) = 0  [extremo realizado M = S]
  d1S <- (b + sigma^2 / 2) * sqrt(T) / sigma
  d2S <- d1S - sigma * sqrt(T)
  
  # Factor potencia (S/K)^{-2b/σ²}: presente cuando d₁ usa log(S/K)
  pot_K <- (S / K)^(-2 * b / sigma^2)
  
  # Factor común del término de corrección
  cf <- S * exp(-r * T) * sigma^2 / (2 * b)
  
  # FIXED CALL:
  
  if (K >= S) {
    # Caso E ≥ M (K ≥ S): d₁ usa log(S/K)
    call <- S * exp(-q * T) * pnorm(d1K) - K * exp(-r * T) * pnorm(d2K) +
      cf * (-pot_K * pnorm(d1K - 2 * b * sqrt(T) / sigma) +
              exp(b * T) * pnorm(d1K))
    
  } else {
    # Caso E < M (K < S): d₁ usa log(S/S) = 0 → d1S
    # Término adicional: (S − K)·e^{-rT} = valor intrínseco descontado
    call <- (S - K) * exp(-r * T) +
      S * exp(-q * T) * pnorm(d1S) - S * exp(-r * T) * pnorm(d2S) +
      cf * (-pnorm(d1S - 2 * b * sqrt(T) / sigma) +
              exp(b * T) * pnorm(d1S))
  }
  
  # FIXED PUT:
  
  if (K <= S) {
    # Caso E ≤ M (K ≤ S): d₁ usa log(S/K)
    put <- K * exp(-r * T) * pnorm(-d2K) - S * exp(-q * T) * pnorm(-d1K) +
      cf * (pot_K * pnorm(-d1K + 2 * b * sqrt(T) / sigma) -
              exp(b * T) * pnorm(-d1K))
    
  } else {
    # Caso E > M (K > S): d₁ usa log(S/S) = 0 → d1S
    # Término adicional: (K − S)·e^{-rT} = valor intrínseco descontado
    put <- (K - S) * exp(-r * T) +
      S * exp(-r * T) * pnorm(-d2S) - S * exp(-q * T) * pnorm(-d1S) +
      cf * (pnorm(-d1S + 2 * b * sqrt(T) / sigma) -
              exp(b * T) * pnorm(-d1S))
  }
  
  list(call = max(call, 0), put = max(put, 0))
}

# --- Barrera ---
# Fórmulas cerradas para opciones barrera estándar.
# Haug, E.G. (2007) Libro "The Complete Guide to Option Pricing Formulas", 2ª ed., sección 4.17.1, p. 152–153.
# Basadas originalmente en Reiner & Rubinstein (1991).

# Válidas para:
#   - Opciones europeas
#   - Monitorización continua de la barrera
#   - Modelo Black-Scholes con dividend yield continuo q
#   - Sin rebate (E = 0, F = 0)

# Convenciones de Haug:
#   phi = +1 (call), -1 (put)
#   eta = +1 (down barriers), -1 (up barriers)
#   b   = r - q  
#   mu  = (b - σ²/2) / σ²

# Las fórmulas tienen dos ramas según la relación entre el strike K
# y la barrera B (notación de Haug: X vs H):
#   - K > B: el strike está por encima de la barrera
#   - K <= B: el strike está por debajo (o igual) a la barrera

# Las validaciones del server garantizan B < S (down) o B > S (up),
# pero K puede estar en cualquier lado respecto a B — de ahí las dos ramas.

bs_barrera <- function(S, K, r, q, sigma, T, B, tipo, tipo_barrera) {
  
  b  <- r - q
  sT <- sigma * sqrt(T)                      # σ√T
  mu <- (b - 0.5 * sigma^2) / sigma^2        # μ = (b - σ²/2) / σ²
  
  phi <- if (tipo == "call") 1L else -1L
  eta <- if (tipo_barrera %in% c("down-and-out", "down-and-in")) 1L else -1L
  
  # Argumentos de la distribución normal (Haug p. 152)
  x1 <- log(S / K) / sT + (1 + mu) * sT
  x2 <- log(S / B) / sT + (1 + mu) * sT
  y1 <- log(B^2 / (S * K)) / sT + (1 + mu) * sT
  y2 <- log(B / S) / sT + (1 + mu) * sT
  
  # Factores de descuento
  Se <- S * exp(-q * T)      # S · e^{-qT}  = S · e^{(b-r)T}
  Ke <- K * exp(-r * T)      # K · e^{-rT}
  
  # Potencias de reflexión (Haug p. 152)
  #   (B/S)^{2(μ+1)}  →  peso del término en S dentro de C y D
  #   (B/S)^{2μ}      →  peso del término en K dentro de C y D
  pow1 <- (B / S)^(2 * (mu + 1))
  pow2 <- (B / S)^(2 * mu)
  
  # Bloques A–D de Haug (E = F = 0 sin rebate)
  A  <- phi * Se * pnorm(phi * x1) - phi * Ke * pnorm(phi * x1 - phi * sT)
  Bb <- phi * Se * pnorm(phi * x2) - phi * Ke * pnorm(phi * x2 - phi * sT)
  C  <- phi * Se * pow1 * pnorm(eta * y1) - phi * Ke * pow2 * pnorm(eta * y1 - eta * sT)
  D  <- phi * Se * pow1 * pnorm(eta * y2) - phi * Ke * pow2 * pnorm(eta * y2 - eta * sT)
  
  # Selección de fórmula según tipo (Haug p. 153)
  # La condición K > B / K <= B es la distinción "X > H" / "X < H" de Haug.
  result <- if (tipo == "call") {
    switch(tipo_barrera,
           "down-and-in"  = if (K > B) C      else A - Bb + D,
           "up-and-in"    = if (K > B) A      else Bb - C + D,
           "down-and-out" = if (K > B) A - C  else Bb - D,
           "up-and-out"   = if (K > B) 0      else A - Bb + C - D,
           NA
    )
  } else {
    switch(tipo_barrera,
           "down-and-in"  = if (K > B) Bb - C + D  else A,
           "up-and-in"    = if (K > B) A - Bb + D  else C,
           "down-and-out" = if (K > B) A - Bb + C - D  else 0,
           "up-and-out"   = if (K > B) Bb - D  else A - C,
           NA
    )
  }
  
  result
}


# -----------------------------------------------------------------------------
# 2.6.3 INTERVALO DE CONFIANZA AL 95% (Teorema Central del Límite)
# -----------------------------------------------------------------------------

# La media muestral de los payoffs Mᵢ sigue aproximadamente una N(µ, σ²/M)
# para M suficientemente grande (TCL). El intervalo para el precio descontado:

#   IC₉₅ = e^{-rT} · [ mean(payoffs) ± z₀.₉₇₅ · sd(payoffs)/√M ]

# El error estándar decrece como 1/√M: para reducirlo a la mitad,
# necesitamos 4 veces más simulaciones. Esta es la "ley de MC".

mc_intervalo_confianza <- function(payoffs, r, T, nivel = 0.95) {
  M    <- length(payoffs)
  mu   <- mean(payoffs)
  s    <- sd(payoffs)              
  z    <- qnorm((1 + nivel) / 2)  
  
  desc <- exp(-r * T)
  
  precio    <- desc * mu
  error_est <- desc * z * s / sqrt(M)
  
  list(
    precio        = precio,
    ic_inf        = precio - error_est,
    ic_sup        = precio + error_est,
    error_estandar = desc * s / sqrt(M),
    M             = M
  )
}


# -----------------------------------------------------------------------------
# 2.7  GRIEGAS MONTE CARLO  (diferencias finitas sobre el precio MC)
# -----------------------------------------------------------------------------
# Para cada griega, se re-ejecutan las simulaciones con parámetro perturbado.
# Se usa la misma semilla → los números aleatorios son idénticos → el error
# estocástico se cancela en la diferencia (técnica de varianza reducida).

# Delta:  Δ ≈ (V(S+h) - V(S-h)) / (2h)           h = 1% de S
# Gamma:  Γ ≈ (V(S+h) - 2V + V(S-h)) / h²
# Vega:   ν ≈ (V(σ+ε) - V(σ-ε)) / (2ε)           ε = 0.001
# Rho:    ρ ≈ (V(r+ε) - V(r-ε)) / (2ε)            ε = 0.0001
# Theta:  Θ ≈ V(T - 1/365) - V(T)                 (cambio en un día)

griegas_mc <- function(S, K, r, q, sigma, T, n_steps, n_sims, tipo_opcion,
                       tipo_mc, semilla,
                       tipo_media=NULL, tipo_barrera=NULL,
                       nivel_barrera=NULL, tipo_strike=NULL) {
  
  sem <- if (!is.null(semilla) && semilla>0) semilla else 111
  M   <- min(n_sims, 3000)   # limitado para velocidad; griega necesita precisión relativa
  
  # Función auxiliar: ejecuta MC y devuelve precio para parámetros dados
  precio_mc <- function(S_v=S, sigma_v=sigma, r_v=r, q_v=q, T_v=T) {
    res <- switch(tipo_mc,
                  "vanilla"  = mc_vanilla(S_v,K,r_v,q_v,sigma_v,T_v,M,tipo_opcion,sem),
                  "asiatica" = mc_asiatica(S_v,K,r_v,q_v,sigma_v,T_v,n_steps,M,tipo_opcion,tipo_media,sem),
                  "barrera"  = mc_barrera(S_v,K,r_v,q_v,sigma_v,T_v,n_steps,M,tipo_opcion,tipo_barrera,nivel_barrera,sem),
                  "lookback" = mc_lookback(S_v,K,r_v,q_v,sigma_v,T_v,n_steps,M,tipo_opcion,tipo_strike,sem)
    )
    exp(-r_v*T_v)*mean(res$payoffs)
  }
  
  h  <- S*0.01;  es <- 0.001;  er <- 0.0001;  eT <- 1/365
  
  V0   <- precio_mc()
  Vsp  <- precio_mc(S_v=S+h);    Vsm  <- precio_mc(S_v=S-h)
  Vsip <- precio_mc(sigma_v=sigma+es); Vsim <- precio_mc(sigma_v=sigma-es)
  Vrp  <- precio_mc(r_v=r+er);    Vrm  <- precio_mc(r_v=r-er)
  
  delta    <- (Vsp-Vsm)/(2*h)
  gamma    <- (Vsp-2*V0+Vsm)/h^2
  vega     <- (Vsip-Vsim)/(2*es)
  rho      <- (Vrp-Vrm)/(2*er)
  
  theta_diario <- NA
  if (T > eT*2) {
    VTm <- precio_mc(T_v=T-eT)
    theta_diario <- VTm - V0   # cambio de valor en un día (típicamente negativo)
  }
  
  list(delta=delta, gamma=gamma, theta_diario=theta_diario,
       vega=vega, vega_pct=vega/100, rho=rho, rho_pct=rho/100)
}


############################################################################### 
# 3. INTERFAZ DE USUARIO (UI)
############################################################################### 

# Aquí definimos el aspecto visual de la aplicación:
# - El panel lateral izquierdo con los inputs del usuario
# - El panel principal derecho con los outputs (resultados y gráficos)

# https://shiny.posit.co/r/getstarted/shiny-basics/lesson1/


ui <- fluidPage(
  
  titlePanel("Interactive Option Pricing Dashboard"),
  
  tags$head(tags$style(HTML("
    body { font-family:'Inter','Segoe UI','Helvetica Neue',Arial,sans-serif; background:#f4f6f8; color:#1f2933; }
    .well { background:#ffffff; border:1px solid #e5e7eb; border-radius:14px; box-shadow:0 4px 14px rgba(15,23,42,0.04); }
    .rbox { background:#ffffff; border:1px solid #e5e7eb; border-radius:14px; padding:18px; margin-bottom:16px; box-shadow:0 6px 18px rgba(15,23,42,0.04); }
    .metric { min-height:118px; }
    .plabel { font-size:13px; color:#697586; margin-bottom:6px; }
    .pbig { font-size:28px; font-weight:700; color:#111827; margin-bottom:4px; }
    .ic { font-size:12px; color:#697586; font-style:italic; }
    .warn { background:#fff7e6; border:1px solid #f0b429; border-radius:12px; padding:14px; margin-bottom:14px; font-size:13px; }
    .info { background:#e8f7fb; border:1px solid #38bdf8; border-radius:12px; padding:14px; margin-bottom:14px; font-size:13px; }
    h4 { color:#1f2933; border-bottom:1px solid #e5e7eb; padding-bottom:8px; margin-top:8px; }
    .nav-tabs > li > a { border-radius:10px 10px 0 0; color:#334e68; }
    .nav-tabs > li.active > a { font-weight:600; border-top:3px solid #2563eb; }
    .btn-primary { background:#2563eb; border-color:#2563eb; border-radius:10px; font-weight:600; }
    .clean-table { width:auto; border-collapse:collapse; font-size:13px; }
    .clean-table th { background:#f3f4f6; font-weight:600; padding:8px 12px; border:1px solid #d1d5db; }
    .clean-table td { padding:7px 12px; border:1px solid #d1d5db; }
    .clean-table tr:nth-child(even) { background:#fafafa; }
  "))),
  
  sidebarLayout(
    sidebarPanel(width = 3,
                 h5("Product", style="font-weight:bold"),
                 selectInput("producto", "Product category",
                             choices = c("Vanilla"="vanilla", "Asian"="asian", "Barrier"="barrier", "Lookback"="lookback"),
                             selected = "vanilla"),
                 selectInput("tipo_opcion", "Option type", choices = c("Call"="call", "Put"="put"), selected = "call"),
                 conditionalPanel("input.producto=='vanilla' || input.producto=='barrier'",
                                  selectInput("ejercicio", "Exercise style", choices = c("European"="european", "American"="american"), selected = "european")
                 ),
                 conditionalPanel("input.producto=='asian' || input.producto=='lookback'",
                                  div(class="info", "Asian and Lookback options are implemented as European-style products in this version.")
                 ),
                 hr(),
                 h5("Market parameters", style="font-weight:bold"),
                 numericInput("S", "Spot price (S0)", 100, min=0.01, step=1),
                 numericInput("K", "Strike (K)", 100, min=0.01, step=1),
                 numericInput("r", "Risk-free rate (r, annual, e.g. 0.05 = 5%)", 0.05, min=0, max=1, step=0.001),
                 
                 # Selector de método de dividendos (solo para opciones vanilla europeas)
                 conditionalPanel(
                   "input.producto == 'vanilla' && input.ejercicio == 'european'",
                   hr(),
                   selectInput("div_metodo", "Dividend method",
                               choices = c(
                                 "Continuous yield (q)" = "continuo",
                                 "Discrete dividends - spot adjustment" = "spot_adj"
                               ),
                               selected = "continuo"
                   ),
                   # Panel para tasa continua q 
                   conditionalPanel(
                     "input.div_metodo == 'continuo'",
                     numericInput("q", "Dividend yield (q, annual, e.g. 0.02 = 2%)", 0, min=0, max=1, step=0.001),
                     p(style="font-size:11px;color:#697586;",
                       "Approximates dividends as a continuous yield. Suitable when there are many small, spread-out dividends.")
                   ),
                   # Panel para dividendos discretos (importes y fechas)
                   conditionalPanel(
                     "input.div_metodo == 'spot_adj'",
                     p(style="font-size:11px;color:#697586;",
                       "Each dividend reduces the initial spot by its present value (spot adjustment method).",
                       br(), "Each time tᵢ must satisfy 0 < tᵢ ≤ T. Each amount Dᵢ must be non-negative."),
                     # Filas de dividendos generadas dinámicamente en el servidor
                     uiOutput("div_inputs_ui"),
                     # Botones para añadir / eliminar la última fila
                     div(
                       class = "form-group",
                       fluidRow(
                         column(6, actionButton("add_div", "+ Add dividend", style="font-size:12px;padding:4px 10px;width:100%;")),
                         column(6, actionButton("remove_div", "- Remove last dividend", style="font-size:12px;padding:4px 10px;width:100%;background:#f3f4f6;color:#374151;border:1px solid #d1d5db;"))
                       )
                     )
                   )
                 ),
                 # Para todos los demás productos input q 
                 conditionalPanel(
                   "!(input.producto == 'vanilla' && input.ejercicio == 'european')",
                   numericInput("q", "Dividend yield (q, annual, e.g. 0.02 = 2%)", 0, min=0, max=1, step=0.001)
                 ),
                 
                 numericInput("sigma", "Volatility (sigma, annual, e.g. 0.20 = 20%)", 0.20, min=0.001, max=5, step=0.01),
                 numericInput("T", "Time to maturity (T, years)", 1, min=0.01, step=0.25),
                 conditionalPanel("input.producto=='vanilla' || input.producto=='barrier'",
                                  hr(),
                                  h5("Binomial parameters", style="font-weight:bold"),
                                  numericInput("N", "Steps (N)", 5, min=2, max=500, step=1),
                                  sliderInput("step_range", "Steps to display in tree", min=0, max=50, value=c(0,5), step=1)
                 ),
                 conditionalPanel("input.producto=='asian'",
                                  hr(),
                                  h5("Asian settings", style="font-weight:bold"),
                                  selectInput("mc_media", "Average type", choices=c("Arithmetic"="aritmetica", "Geometric"="geometrica"), selected="aritmetica")
                 ),
                 conditionalPanel("input.producto=='barrier'",
                                  hr(),
                                  h5("Barrier settings", style="font-weight:bold"),
                                  selectInput("mc_tipo_barrera", "Barrier type",
                                              choices=c("Down-and-Out"="down-and-out", "Down-and-In"="down-and-in", "Up-and-Out"="up-and-out", "Up-and-In"="up-and-in"),
                                              selected="down-and-out"),
                                  numericInput("mc_nivel_barrera", "Barrier level (B)", 90, min=0.01, step=1),
                                  p(style="font-size:11px;color:#697586;", "Down barriers are usually below S0. Up barriers are usually above S0.")
                 ),
                 conditionalPanel("input.producto=='lookback'",
                                  hr(),
                                  h5("Lookback settings", style="font-weight:bold"),
                                  selectInput("mc_lookback_strike", "Strike type",
                                              choices=c("Floating"="floating", "Fixed"="fixed"), selected="floating")
                 ),
                 br(),
                 actionButton("calcular", "Calculate", class="btn-primary btn-block")
    ),
    mainPanel(width=9, uiOutput("tabs_main"))
  )
)


############################################################################### 
# 4. SERVIDOR
############################################################################### 

# El servidor contiene la lógica reactiva de Shiny:
# - Lee los inputs del usuario
# - Ejecuta las funciones matemáticas
# - Renderiza los outputs (tablas, gráficos, textos)

# Los outputs se actualizan automáticamente
# cada vez que cambia un input que usan. Aquí usamos eventReactive()


server <- function(input, output, session) {
  
  app_started <- reactiveVal(FALSE)
  
  mc_params <- reactiveVal(list(
    n_sims  = 10000,
    n_steps = 252,
    semilla = 111
  ))
  
  observeEvent(input$calcular, {
    app_started(TRUE)
    
    old <- mc_params()
    mc_params(list(
      n_sims  = input$mc_n_sims  %||% old$n_sims,
      n_steps = input$mc_n_steps %||% old$n_steps,
      semilla = input$mc_semilla %||% old$semilla
    ))
  })
  
  observeEvent(input$update_mc, {
    old <- mc_params()
    mc_params(list(
      n_sims  = input$mc_n_sims  %||% old$n_sims,
      n_steps = input$mc_n_steps %||% old$n_steps,
      semilla = input$mc_semilla %||% old$semilla
    ))
  })
  
  # Gestión dinámica de filas de dividendos discretos
  MAX_DIVS <- 10
  n_divs <- reactiveVal(1)  # por defecto 1 fila visible
  
  # Permite añadir una fila si no se han añadido ya las 10
  observeEvent(input$add_div, {
    if (n_divs() < MAX_DIVS) n_divs(n_divs() + 1)
  })
  
  # Eliminar la última fila si hay más de una
  observeEvent(input$remove_div, {
    if (n_divs() > 1) n_divs(n_divs() - 1)
  })
  
  output$div_inputs_ui <- renderUI({
    n <- n_divs()
    rows <- lapply(seq_len(n), function(i) {
      # Creamos el nombre que se verá para cada dividendo: D₁, D₂, D₃...
      # Si el número es del 1 al 9, usamos el formato con subíndice.
      # Si llega a 10, lo escribimos directamente como "10" para evitar problemas.
      sub_char <- if (i <= 9) {
        intToUtf8(0x2080 + i) 
      } else {
        "10"
      }
      fluidRow(
        column(6,
               numericInput(paste0("div", i, "_imp"),
                            paste0("D", sub_char, " (amount)"),
                            value = isolate(input[[paste0("div", i, "_imp")]]) %||% 0,
                            min = 0, step = 0.5)
        ),
        column(6,
               numericInput(paste0("div", i, "_t"),
                            paste0("t", sub_char, " (years)"),
                            value = isolate(input[[paste0("div", i, "_t")]]) %||% round(i * 0.5, 1),
                            min = 0.001, step = 0.1)
        )
      )
    })
    
    # Indicador del límite máximo cuando se ha alcanzado
    limite_msg <- if (n >= MAX_DIVS) {
      p(style="font-size:11px;color:#f97316;margin-top:4px;",
        paste0("Maximum of ", MAX_DIVS, " dividends reached."))
    } else NULL
    
    tagList(rows, limite_msg)
  })
  
  # El rango máximo del slider se adapta al número de pasos elegido
  observe({
    max_steps <- min(input$N %||% 50, 50)
    current <- input$step_range %||% c(0, min(5, max_steps))
    current[1] <- max(0, min(current[1], max_steps))
    current[2] <- max(current[1], min(current[2], max_steps))
    updateSliderInput(session, "step_range", min=0, max=max_steps, value=current)
  })
  
  ej_ef <- reactive({
    if (input$producto %in% c("asian", "lookback")) "european" else input$ejercicio
  })
  
  mc_inputs <- reactive({
    mc_params()
  })
  
  clean_table <- function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    tags$table(class="clean-table",
               tags$thead(tags$tr(lapply(names(df), tags$th))),
               tags$tbody(lapply(seq_len(nrow(df)), function(i) {
                 tags$tr(lapply(df[i, ], function(x) tags$td(as.character(x))))
               }))
    )
  }
  
  make_edges <- function(n_min, n_max) {
    edges <- data.frame()
    if (n_max <= n_min) return(edges)
    for (n in n_min:(n_max - 1)) {
      for (j in 0:n) {
        x0 <- n; y0 <- j - n / 2
        edges <- rbind(edges,
                       data.frame(x=x0, y=y0, xend=n+1, yend=(j+1)-(n+1)/2, move="Up"),
                       data.frame(x=x0, y=y0, xend=n+1, yend=j-(n+1)/2, move="Down")
        )
      }
    }
    edges
  }
  
  payoff_plot <- function(S, K, tipo) {
    x_min <- max(0.01, min(S, K) * 0.4)
    x_max <- max(S, K) * 1.8
    x <- seq(x_min, x_max, length.out=250)
    payoff <- if (tipo == "call") pmax(x - K, 0) else pmax(K - x, 0)
    plot_ly(x=x, y=payoff, type="scatter", mode="lines",
            line=list(color="#2563eb", width=3),
            hovertemplate="Underlying price: %{x:.4f}<br>Payoff: %{y:.4f}<extra></extra>") %>%
      layout(title="Payoff at maturity",
             xaxis=list(title="Underlying price at maturity"),
             yaxis=list(title="Payoff"),
             shapes=list(list(type="line", x0=K, x1=K, y0=0, y1=max(payoff), line=list(color="#ef4444", dash="dash"))),
             margin=list(l=50,r=20,b=50,t=50))
  }
  
  rv <- reactive({
    req(app_started())
    
    # Para opciones vanilla europeas, el usuario puede elegir entre la tasa continua q 
    # o el ajuste del spot inicial con dividendos discretos 
    # Para cualquier otro producto, se usa siempre la tasa continua q
    
    div_metodo_activo <- if (!is.null(input$div_metodo) &&
                             input$producto == "vanilla" &&
                             (input$ejercicio %||% "european") == "european") {
      input$div_metodo
    } else {
      "continuo"
    }
    
    # Si se usa tasa continua, leemos el input para q; si es spot_adj, q=0
    q_efectivo <- if (div_metodo_activo == "spot_adj") {
      0
    } else {
      input$q %||% 0
    }
    
    # Lectura dinámica de dividendos:
    T_val <- input$T %||% 1   # necesitamos T para la validación de tiempos
    
    divs_discretos <- if (div_metodo_activo == "spot_adj") {
      n <- n_divs()
      importes <- sapply(seq_len(n), function(i) input[[paste0("div", i, "_imp")]] %||% 0)
      tiempos  <- sapply(seq_len(n), function(i) input[[paste0("div", i, "_t")]]  %||% (i * 0.5))
      data.frame(importe = importes, tiempo = tiempos)
    } else {
      data.frame(importe = numeric(0), tiempo = numeric(0))
    }
    
    validate(
      need(input$S > 0, "Spot must be positive"),
      need(input$K > 0, "Strike must be positive"),
      need(q_efectivo >= 0, "Dividend yield must be zero or positive"),
      need(input$sigma > 0, "Volatility must be positive"),
      need(T_val > 0, "Time to maturity must be positive")
    )
    
    if (div_metodo_activo == "spot_adj" && nrow(divs_discretos) > 0) {
      
      # Comprobamos que los importes no sean negativos
      neg_idx <- which(divs_discretos$importe < 0)
      validate(
        need(
          length(neg_idx) == 0,
          paste0("The amount of dividend D", neg_idx[1], " is negative (",
                 divs_discretos$importe[neg_idx[1]], "). ",
                 "Dividends must be positive amounts or zero.")
        )
      )
      
      # Comprobamos que todos los tiempos estén en el intervalo (0, T]
      fuera_idx <- which(divs_discretos$tiempo <= 0 | divs_discretos$tiempo > T_val)
      validate(
        need(
          length(fuera_idx) == 0,
          paste0("The time t", fuera_idx[1], " = ", divs_discretos$tiempo[fuera_idx[1]],
                 " years is outside the valid interval (0, T] = (0, ", T_val, "]. ",
                 "Enter a strictly positive time less than or equal to the maturity T.")
        )
      )
      
      # Comprobamos que el valor presente de los dividendos no supere el precio inicial
      pv_divs <- sum(divs_discretos$importe * exp(-input$r * divs_discretos$tiempo))
      
      validate(
        need(
          pv_divs < input$S,
          paste0("The present value of the dividends is greater than or equal to the spot price. ",
                 "This would make the adjusted spot S* non-positive. ",
                 "Please reduce the dividend amounts or check the input values.")
        )
      )
    }
    
    S <- input$S; K <- input$K; r <- input$r; q <- q_efectivo
    sigma <- input$sigma; T <- input$T
    tipo <- input$tipo_opcion; prod <- input$producto; ej <- ej_ef()
    out <- list(prod=prod, ej=ej, tipo=tipo, S=S, K=K, r=r, q=q, sigma=sigma, T=T,
                div_metodo=div_metodo_activo, divs_discretos=divs_discretos)
    mci_actual <- mc_inputs()
    out$mc_n_sims  <- mci_actual$n_sims
    out$mc_n_steps <- mci_actual$n_steps
    out$mc_semilla <- mci_actual$semilla
    
    if (prod %in% c("vanilla", "barrier")) {
      validate(need(input$N >= 2, "N must be at least 2"))
      out$N <- input$N
      sr <- input$step_range %||% c(0, min(5, input$N))
      out$n_min <- max(0, min(sr[1], input$N))
      out$n_max <- max(out$n_min, min(sr[2], input$N))
      if (out$n_max < out$n_min) out$n_max <- out$n_min
      
      # La condición necesaria para que p* ∈ [0,1] es d < exp((r-q)·δt) < u
      # Si no se cumple, el árbol no es válido y los resultados carecen de sentido.
      dt_check <- T / input$N
      u_check  <- exp(sigma * sqrt(dt_check))
      d_check  <- 1 / u_check
      p_check  <- (exp((r - q) * dt_check) - d_check) / (u_check - d_check)
      out$p_check <- p_check  
      
      validate(
        need(
          is.finite(p_check) && p_check >= 0 && p_check <= 1,
          paste0(
            "\u26a0 The risk-neutral probability is outside the interval [0,1] ",
            "(p* = ", round(p_check, 4), "), so the binomial tree is not valid ",
            "with these parameters. Try increasing the number of tree steps ",
            "to reduce the interest rate per period (r*deltat), or check the values ",
            "entered for 'r' and 'q'."
          )
        )
      )
    }
    
    if (prod == "vanilla") {
      
      # Bifurcación según el método de dividendos elegido
      
      if (div_metodo_activo == "spot_adj") {
        
        res_bin <- binomial_model_spot_adj(S, K, r, sigma, T, out$N, tipo, ej, divs_discretos)
        out$bin <- res_bin
        out$arbol <- preparar_datos_arbol(res_bin, K, tipo, ej)
        out$gr_bin <- griegas_binomiales(res_bin, res_bin$S_star, K, r, q=0, sigma, T, out$N, tipo, ej)
        
        if (ej == "european") {
          # Black-Scholes con S* ajustado
          out$bs <- black_scholes_spot_adj(S, K, r, sigma, T, tipo, divs_discretos)
          out$gr_bs <- griegas_bs(res_bin$S_star, K, r, q=0, sigma, T, tipo)
          out$S_star <- res_bin$S_star  # guardamos S* para mostrarlo en la interfaz
          
          mci <- mc_inputs(); nsims <- mci$n_sims; nst <- mci$n_steps; sem <- mci$semilla
          validate(need(nsims >= 100, "At least 100 Monte Carlo simulations are required"))
          
          # Monte Carlo con S* ajustado y q=0
          res_mc <- mc_vanilla_spot_adj(S, K, r, sigma, T, nsims, tipo, sem, divs_discretos)
          out$mc <- res_mc; out$ic <- mc_intervalo_confianza(res_mc$payoffs, r, T)
          out$mc_ref <- list(precio=out$bs$precio, etiqueta="Black-Scholes with spot adjustment (analytical reference)")
          out$gr_mc  <- griegas_mc(res_bin$S_star, K, r, q=0, sigma, T, 1, nsims, tipo, "vanilla", sem)
          
          # Tabla GBM vs Euler (usando S*)
          nc <- min(nsims, 3000); s2 <- if(sem>0) sem else 111
          rex <- mc_simular_gbm(res_bin$S_star, r, q=0, sigma, T, nst, nc, s2)
          pex <- if(tipo=="call") pmax(rex$caminos[,nst+1]-K,0) else pmax(K-rex$caminos[,nst+1],0)
          dt_e <- T/nst; ce <- matrix(0, nc, nst+1); ce[,1] <- res_bin$S_star
          for(ps in 1:nst) ce[,ps+1] <- ce[,ps] + r*ce[,ps]*dt_e + sigma*ce[,ps]*sqrt(dt_e)*rex$phi[,ps]
          peu <- if(tipo=="call") pmax(ce[,nst+1]-K,0) else pmax(K-ce[,nst+1],0)
          r1 <- mc_simular_gbm(res_bin$S_star, r, q=0, sigma, T, 1, nc, s2)
          p1 <- if(tipo=="call") pmax(r1$caminos[,2]-K,0) else pmax(K-r1$caminos[,2],0)
          out$euler_df <- data.frame(
            Method=c(paste0("GBM exact, 1 step (spot adj., S*=",round(res_bin$S_star,4),")"),
                     paste0("GBM exact, ",nst," steps (spot adj.)"),
                     paste0("Euler, ",nst," steps (spot adj.)"),
                     "Black-Scholes with spot adjustment"),
            Price=round(c(exp(-r*T)*mean(p1), exp(-r*T)*mean(pex),
                          exp(-r*T)*mean(peu), out$bs$precio), 4),
            Note=c("No discretisation error", "No discretisation error",
                   "Euler discretisation error", "Analytical solution")
          )
        }
        
      } else {
        
        # Método continuo (tasa q)
        res_bin <- binomial_model(S,K,r,q,sigma,T,out$N,tipo,ej)
        out$bin <- res_bin
        out$arbol <- preparar_datos_arbol(res_bin, K, tipo, ej)
        out$gr_bin <- griegas_binomiales(res_bin,S,K,r,q,sigma,T,out$N,tipo,ej)
        
        if (ej == "european") {
          out$bs <- black_scholes(S,K,r,q,sigma,T,tipo)
          out$gr_bs <- griegas_bs(S,K,r,q,sigma,T,tipo)
          mci <- mc_inputs(); nsims <- mci$n_sims; nst <- mci$n_steps; sem <- mci$semilla
          validate(need(nsims >= 100, "At least 100 Monte Carlo simulations are required"))
          res_mc <- mc_vanilla(S,K,r,q,sigma,T,nsims,tipo,sem)
          out$mc <- res_mc; out$ic <- mc_intervalo_confianza(res_mc$payoffs,r,T)
          out$mc_ref <- list(precio=out$bs$precio, etiqueta="Black-Scholes analytical reference")
          out$gr_mc <- griegas_mc(S,K,r,q,sigma,T,1,nsims,tipo,"vanilla",sem)
          
          nc <- min(nsims,3000); s2 <- if(sem>0) sem else 111
          rex <- mc_simular_gbm(S,r,q,sigma,T,nst,nc,s2)
          pex <- if(tipo=="call") pmax(rex$caminos[,nst+1]-K,0) else pmax(K-rex$caminos[,nst+1],0)
          dt_e <- T/nst; ce <- matrix(0,nc,nst+1); ce[,1] <- S
          for(ps in 1:nst) ce[,ps+1] <- ce[,ps] + (r - q)*ce[,ps]*dt_e + sigma*ce[,ps]*sqrt(dt_e)*rex$phi[,ps]
          peu <- if(tipo=="call") pmax(ce[,nst+1]-K,0) else pmax(K-ce[,nst+1],0)
          r1 <- mc_simular_gbm(S,r,q,sigma,T,1,nc,s2)
          p1 <- if(tipo=="call") pmax(r1$caminos[,2]-K,0) else pmax(K-r1$caminos[,2],0)
          out$euler_df <- data.frame(
            Method=c(paste0("GBM exact, 1 step, dt = T = ",T), paste0("GBM exact, ",nst," steps"), paste0("Euler, ",nst," steps, same shocks"), "Black-Scholes reference"),
            Price=round(c(exp(-r*T)*mean(p1),exp(-r*T)*mean(pex),exp(-r*T)*mean(peu),out$bs$precio),4),
            Note=c("No discretisation error", "No discretisation error", "Euler discretisation error", "Analytical solution")
          )
        }
      }  
      
    } else if (prod == "barrier") {
      B <- input$mc_nivel_barrera
      tb <- input$mc_tipo_barrera
      
      validate(
        need(B > 0, "Barrier level must be positive"),
        need(B != S, "Barrier level must differ from spot"),
        need(
          !(tb %in% c("down-and-out", "down-and-in") && B >= S),
          "For a down barrier, the barrier level must be below the spot price."
        ),
        need(
          !(tb %in% c("up-and-out", "up-and-in") && B <= S),
          "For an up barrier, the barrier level must be above the spot price."
        )
      )
      out$B <- B
      out$tb <- tb
      res_bar <- binomial_barrier_model(S,K,r,q,sigma,T,out$N,tipo,ej,tb,B)
      out$bar <- res_bar
      out$arbol_bar <- preparar_datos_arbol_barrera(res_bar, K, tipo, ej)
      out$gr_bar <- griegas_binomiales_barrera(res_bar,S,K,r,q,sigma,T,out$N,tipo,ej,tb,B)
      if (ej == "european") {
        mci <- mc_inputs(); nsims <- mci$n_sims; nst <- mci$n_steps; sem <- mci$semilla
        validate(need(nsims >= 100, "At least 100 Monte Carlo simulations are required"))
        res_mc <- mc_barrera(S,K,r,q,sigma,T,nst,nsims,tipo,tb,B,sem)
        out$mc <- res_mc; out$ic <- mc_intervalo_confianza(res_mc$payoffs,r,T)
        vbs <- bs_barrera(S,K,r,q,sigma,T,B,tipo,tb)
        out$mc_ref <- list(precio=vbs, etiqueta="Analytical barrier reference, continuous monitoring")
        out$gr_mc <- griegas_mc(S,K,r,q,sigma,T,nst,nsims,tipo,"barrera",sem,tipo_barrera=tb,nivel_barrera=B)
      }
    } else if (prod == "asian") {
      mci <- mc_inputs(); nsims <- mci$n_sims; nst <- mci$n_steps; sem <- mci$semilla
      tm <- input$mc_media %||% "aritmetica"
      validate(need(nsims >= 100, "At least 100 Monte Carlo simulations are required"))
      res_mc <- mc_asiatica(S,K,r,q,sigma,T,nst,nsims,tipo,tm,sem)
      out$mc <- res_mc; out$ic <- mc_intervalo_confianza(res_mc$payoffs,r,T); out$tipo_media <- tm
      if (tm == "geometrica") {
        ag <- bs_asiatica_geometrica(S,K,r,q,sigma,T)
        out$mc_ref <- list(precio=if(tipo=="call") ag$call else ag$put, etiqueta="Geometric Asian analytical reference")
      } else {
        ar <- bs_asiatica_aritmetica(S,K,r,q,sigma,T)
        out$mc_ref <- list(precio=if(tipo=="call") ar$call else ar$put, etiqueta="Arithmetic Asian moment-matching approximation")
      }
      out$gr_mc <- griegas_mc(S,K,r,q,sigma,T,nst,nsims,tipo,"asiatica",sem,tipo_media=tm)
    } else if (prod == "lookback") {
      mci <- mc_inputs(); nsims <- mci$n_sims; nst <- mci$n_steps; sem <- mci$semilla
      ts <- input$mc_lookback_strike %||% "floating"
      validate(need(nsims >= 100, "At least 100 Monte Carlo simulations are required"))
      res_mc <- mc_lookback(S,K,r,q,sigma,T,nst,nsims,tipo,ts,sem)
      out$mc <- res_mc; out$ic <- mc_intervalo_confianza(res_mc$payoffs,r,T); out$tipo_strike <- ts
      if (ts == "floating") {
        lb <- bs_lookback_floating(S,r,q,sigma,T)
        out$mc_ref <- list(precio=if(tipo=="call") lb$call else lb$put, etiqueta="Lookback floating analytical reference")
      } else {
        lb <- bs_lookback_fixed(S,K,r,q,sigma,T)
        out$mc_ref <- list(precio=if(tipo=="call") lb$call else lb$put, etiqueta="Lookback fixed analytical reference")
      }
      out$gr_mc <- griegas_mc(S,K,r,q,sigma,T,nst,nsims,tipo,"lookback",sem,tipo_strike=ts)
    }
    out
  })
  
  output$tabs_main <- renderUI({
    if (!app_started()) {
      return(div(class="rbox", style="padding:50px;text-align:center;color:#697586;",
                 h4("Configure the contract and press Calculate"),
                 p("After the first calculation, changes in the inputs update the results automatically.")))
    }
    req(rv())
    r0 <- rv(); p0 <- r0$prod; e0 <- r0$ej
    mk <- function(...) {tabsetPanel(id = "tabs_inner", selected = isolate(input$tabs_inner), ...)}
    if (p0=="vanilla" && e0=="european") {
      mk(tabPanel("Overview", br(), uiOutput("tab_overview")), tabPanel("Binomial Tree", br(), uiOutput("tab_tree")), tabPanel("Greeks", br(), uiOutput("tab_greeks")), tabPanel("Monte Carlo", br(), uiOutput("tab_mc")), tabPanel("MC Greeks", br(), uiOutput("tab_mc_greeks")))
    } else if (p0=="vanilla" && e0=="american") {
      mk(tabPanel("Overview", br(), uiOutput("tab_overview")), tabPanel("Binomial Tree", br(), uiOutput("tab_tree")), tabPanel("Greeks", br(), uiOutput("tab_greeks")))
    } else if (p0=="barrier" && e0=="european") {
      mk(tabPanel("Overview", br(), uiOutput("tab_overview")), tabPanel("Barrier Tree", br(), uiOutput("tab_barrier_tree")), tabPanel("Greeks", br(), uiOutput("tab_greeks")), tabPanel("Monte Carlo", br(), uiOutput("tab_mc")), tabPanel("MC Greeks", br(), uiOutput("tab_mc_greeks")))
    } else if (p0=="barrier" && e0=="american") {
      mk(tabPanel("Overview", br(), uiOutput("tab_overview")), tabPanel("Barrier Tree", br(), uiOutput("tab_barrier_tree")), tabPanel("Greeks", br(), uiOutput("tab_greeks")))
    } else {
      mk(tabPanel("Monte Carlo", br(), uiOutput("tab_mc")), tabPanel("MC Greeks", br(), uiOutput("tab_mc_greeks")))
    }
  })
  
  output$tab_overview <- renderUI({
    req(rv()); r0 <- rv(); p0 <- r0$prod; e0 <- r0$ej
    if (p0 == "vanilla") {
      # Mostramos p* en la tabla de parámetros binomiales
      params_df <- data.frame(
        Parameter=c("N","dt","u","d","p*"),
        Formula=c("N","T/N","exp(sigma sqrt(dt))","1/u","(exp((r-q) dt)-d)/(u-d)"),
        Value=round(c(r0$N,r0$bin$dt,r0$bin$u,r0$bin$d,r0$bin$p),6)
      )
      
      div_info_panel <- if (!is.null(r0$div_metodo) && r0$div_metodo == "spot_adj" && !is.null(r0$divs_discretos) && sum(r0$divs_discretos$importe, na.rm = TRUE) > 0) {
        S_star_val <- if (!is.null(r0$S_star)) round(r0$S_star, 4) else round(r0$bin$S_star, 4)
        div(class="info",
            strong("Dividend method: Discrete dividends - spot adjustment"),
            br(),
            paste0("Adjusted spot S* = S₀ - PV(dividends) = ", r0$S, " - ",
                   round(r0$S - S_star_val, 4), " = ", S_star_val),
            br(),
            "All three pricing methods (Binomial, Black-Scholes, Monte Carlo) use S* with q = 0.",
            br(),
            em("Note: this approximation is only reliable for European vanilla options (not American, not path-dependent).")
        )
      } else {
        NULL
      }
      
      tagList(
        fluidRow(
          column(4, div(class="rbox metric", p(class="plabel","Binomial price"), p(class="pbig", round(r0$bin$precio,4)))),
          column(4, div(class="rbox metric", p(class="plabel","Black-Scholes"), if(e0=="european") tagList(p(class="pbig", round(r0$bs$precio,4))) else div(class="warn","Black-Scholes is not valid for American options."))),
          column(4, div(class="rbox metric", p(class="plabel","Difference"), if(e0=="european") { dif<-r0$bin$precio-r0$bs$precio; p(class="pbig", style=paste0("color:", if(abs(dif)<0.01) "#16a34a" else "#f97316"), round(dif,4)) } else p(style="color:#697586","N/A")))
        ),
        div_info_panel,
        h4("Binomial parameters"), div(class="rbox", clean_table(params_df)),
        h4("Payoff diagram"), div(class="rbox", plotlyOutput("payoff_plot", height="320px"))
      )
    } else if (p0 == "barrier") {
      params_df <- data.frame(
        Parameter=c("N","dt","u","d","p*"),
        Value=c(round(c(r0$N,r0$bar$dt,r0$bar$u,r0$bar$d,r0$bar$p),6))
      )
      
      v_ref <- if (e0 == "european") {bs_barrera(r0$S, r0$K, r0$r, r0$q, r0$sigma, r0$T, r0$B, r0$tipo, r0$tb)
      } else {NA}
      
      dif_bar <- r0$bar$precio - v_ref
      
      tagList(
        fluidRow(
          column(4, div(class="rbox metric", p(class="plabel",paste("Binomial price")), p(class="pbig", round(r0$bar$precio,4)))),
          column(4, div(class="rbox metric", p(class="plabel","Analytical reference"), if(e0=="european") { v<-bs_barrera(r0$S,r0$K,r0$r,r0$q,r0$sigma,r0$T,r0$B,r0$tipo,r0$tb); if(!is.na(v)) tagList(p(class="pbig", round(v,4)), p(class="ic","European, continuous monitoring")) else p(style="color:#697586","Not available") } else div(class="warn","Analytical barrier formulas are not used for American barrier options."))),
          column(4, div(class="rbox metric", p(class = "plabel", "Difference"), if (e0 == "european" && !is.na(v_ref)) { tagList(
            p(class = "pbig", style = paste0("color:", if (abs(dif_bar) < 0.10) "#16a34a" else "#f97316"), round(dif_bar, 4)), p(class = "ic", "Binomial - analytical"))
          } else {
            tagList(p(class = "pbig", style = "color:#697586;", "N/A"), p(class = "ic", "No analytical comparison"))}))),
        h4("Binomial barrier parameters"), div(class="rbox", clean_table(params_df))
      )
    }
  })
  
  output$payoff_plot <- renderPlotly({
    req(rv()); r0 <- rv()
    validate(need(r0$prod == "vanilla", "Payoff plot is shown for vanilla options only."))
    payoff_plot(r0$S, r0$K, r0$tipo)
  })
  
  output$tab_tree <- renderUI({
    req(rv(), rv()$arbol)
    tagList(
      h4("Binomial Tree"),
      p(style="color:#697586;font-size:13px;", if (rv()$ej == "american") "Hover over nodes to see the underlying price, option value, intrinsic value, continuation value and early-exercise information." else "Hover over nodes to see the underlying price and option value."),
      plotlyOutput("arbol_plot", height="560px")
    )
  })
  
  output$arbol_plot <- renderPlotly({
    req(rv(), rv()$arbol)
    r0 <- rv(); datos <- r0$arbol
    n_min <- r0$n_min; n_max <- r0$n_max
    datos <- datos[datos$paso >= n_min & datos$paso <= n_max, ]
    edges <- make_edges(n_min, n_max)
    
    p <- plot_ly()
    if (nrow(edges) > 0) {
      for (i in seq_len(nrow(edges))) {
        p <- add_segments(p, x=edges$x[i], xend=edges$xend[i], y=edges$y[i], yend=edges$yend[i],
                          line=list(color="#9ca3af", width=1.2), hoverinfo="none", showlegend=FALSE)
      }
    }
    
    datos$grupo <- ifelse(datos$ITM, "ITM", "OTM")
    color_map <- c("ITM" = "#008B45", "OTM" = "#8B3626")
    datos_normales <- datos[!datos$early_exercise, ]
    datos_early <- datos[datos$early_exercise, ]
    
    # Nodos normales, sin borde especial
    if (nrow(datos_normales) > 0) {
      p <- add_markers(p, data = datos_normales, x = ~paso, y = ~y_pos, color = ~grupo, colors = color_map,
                       marker = list(size = 18, line = list(width = 0, color = "rgba(0,0,0,0)")), text = ~hover, hoverinfo = "text", showlegend = TRUE)}
    
    # Nodos con ejercicio anticipado, misma lógica de color pero borde rojo
    if (nrow(datos_early) > 0) {
      p <- add_markers(p, data = datos_early,x = ~paso, y = ~y_pos, color = ~grupo, colors = color_map,
                       marker = list(size = 22, line = list(width = 4, color = "#dc2626")), text = ~hover, hoverinfo = "text", showlegend = FALSE)}
    p <- add_text(p, data=datos, x=~paso, y=~(y_pos+0.22), text=~round(S,2), textfont=list(size=11, color="#111827"), showlegend=FALSE, hoverinfo="none")
    layout(p, xaxis=list(title="Time step", dtick=1, zeroline = FALSE, showline = FALSE, range = list(n_min - 0.5, n_max + 0.5)), yaxis=list(title = "", showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, ticks = "", visible = FALSE),
           legend=list(orientation="v"), margin=list(l=40,r=20,b=50,t=30)) %>%
      config(displaylogo=FALSE, modeBarButtonsToRemove=c("select2d","lasso2d","autoScale2d","hoverClosestCartesian","hoverCompareCartesian"))
  })
  
  output$tab_barrier_tree <- renderUI({
    req(rv(), rv()$arbol_bar)
    tagList(
      h4(paste("Barrier Tree", toupper(rv()$tb), "B =", rv()$B)),
      p(style="color:#697586;font-size:13px;", if (rv()$ej == "american") "Hover over nodes to see price, option value, intrinsic value, continuation value, barrier state and early-exercise information." else "Hover over nodes to see price, option value and barrier state."),
      plotlyOutput("arbol_bar_plot", height="560px")
    )
  })
  
  output$arbol_bar_plot <- renderPlotly({
    req(rv(), rv()$arbol_bar)
    r0 <- rv(); datos <- r0$arbol_bar
    n_min <- r0$n_min; n_max <- r0$n_max
    datos <- datos[datos$paso >= n_min & datos$paso <= n_max, ]
    edges <- make_edges(n_min, n_max)
    
    p <- plot_ly()
    if (nrow(edges) > 0) {
      for (i in seq_len(nrow(edges))) {
        p <- add_segments(p, x = edges$x[i], xend = edges$xend[i],
                          y = edges$y[i], yend = edges$yend[i],
                          line = list(color = "#9ca3af", width = 1.2),
                          hoverinfo = "none", showlegend = FALSE)
      }
    }
    
    color_map <- c("Alive" = "#008B45", "Knocked out" = "#8B3626",
                   "Activated" = "#008B45", "Not activated" = "#6b7280")
    
    datos_normales <- datos[!datos$early_exercise, ]
    datos_early <- datos[datos$early_exercise, ]
    
    if (nrow(datos_normales) > 0) {
      p <- add_markers(p, data = datos_normales, x = ~paso, y = ~y_pos,
                       color = ~estado, colors = color_map,
                       marker = list(size = 18, line = list(width = 0, color = "rgba(0,0,0,0)")),
                       text = ~hover, hoverinfo = "text", showlegend = TRUE)
    }
    
    if (nrow(datos_early) > 0) {
      p <- add_markers(p, data = datos_early, x = ~paso, y = ~y_pos,
                       color = ~estado, colors = color_map,
                       marker = list(size = 23, line = list(width = 4, color = "#dc2626")),
                       text = ~hover, hoverinfo = "text", showlegend = FALSE)
    }
    
    p <- add_text(p, data = datos, x = ~paso, y = ~(y_pos + 0.22),
                  text = ~round(S, 2), textfont = list(size = 11, color = "#111827"),
                  showlegend = FALSE, hoverinfo = "none")
    
    layout(p,
           xaxis = list(title = "Time step", dtick = 1, zeroline = FALSE, showline = FALSE, range = list(n_min - 0.5, n_max + 0.5)),
           yaxis = list(title = "", showticklabels = FALSE, showgrid = FALSE,
                        zeroline = FALSE, showline = FALSE, ticks = "", visible = FALSE),
           legend = list(orientation = "v"),
           margin = list(l = 40, r = 20, b = 50, t = 30)) %>%
      config(displaylogo = FALSE,
             modeBarButtonsToRemove = c("select2d", "lasso2d", "autoScale2d",
                                        "hoverClosestCartesian", "hoverCompareCartesian"))
  })
  
  output$tab_greeks <- renderUI({
    req(rv()); r0 <- rv(); p0 <- r0$prod; e0 <- r0$ej
    greek_box <- function(label, value, formula, desc="") {
      div(class="rbox", h5(label), p(style="font-size:20px;font-weight:bold;", ifelse(is.na(value), "N/A", round(as.numeric(value),6))), p(style="font-size:11px;color:#697586;", formula), if(nchar(desc)>0) p(style="font-size:11px;color:#697586;", desc))
    }
    gr_bin <- if (p0=="vanilla") r0$gr_bin else r0$gr_bar
    bin_panel <- tagList(
      h4(if(p0=="barrier") "Binomial barrier Greeks" else "Binomial Greeks"),
      fluidRow(
        column(3, greek_box("Delta", gr_bin$delta, "First sensitivity to S")),
        column(3, greek_box("Gamma", gr_bin$gamma, "Sensitivity of Delta to S")),
        column(3, greek_box("Theta", gr_bin$theta_diario, "Approximate daily time decay")),
        column(3, greek_box("Vega", gr_bin$vega_pct, "Sensitivity to a 1% volatility change"))
      ),
      fluidRow(column(3, greek_box("Rho", gr_bin$rho_pct, "Sensitivity to a 1% rate change")))
    )
    if (p0=="vanilla" && e0=="european") {
      gr_bs <- r0$gr_bs
      df <- data.frame(
        Greek=c("Delta","Gamma","Theta/day","Vega/1% sigma","Rho/1% r"),
        Binomial=round(c(gr_bin$delta,gr_bin$gamma,gr_bin$theta_diario,gr_bin$vega_pct,gr_bin$rho_pct),6),
        Black_Scholes=round(c(gr_bs$delta,gr_bs$gamma,gr_bs$theta_diario,gr_bs$vega_pct,gr_bs$rho_pct),6)
      )
      df$Difference <- round(df$Binomial - df$Black_Scholes, 6)
      tagList(
        h4("Black-Scholes Greeks"),
        fluidRow(
          column(3, greek_box("Delta", gr_bs$delta, "Analytical Delta")),
          column(3, greek_box("Gamma", gr_bs$gamma, "Analytical Gamma")),
          column(3, greek_box("Theta/day", gr_bs$theta_diario, "Analytical daily Theta")),
          column(3, greek_box("Vega", gr_bs$vega_pct, "Analytical Vega per 1%"))
        ),
        fluidRow(column(3, greek_box("Rho", gr_bs$rho_pct, "Analytical Rho per 1%"))),
        bin_panel,
        h4("Comparison"), div(class="rbox", clean_table(df))
      )
    } else {
      tagList(div(class="info", if(p0=="barrier") "Black-Scholes vanilla Greeks are not shown for barrier options." else "Black-Scholes Greeks are not applicable to American options."), bin_panel)
    }
  })
  
  
  output$tab_mc <- renderUI({req(rv()); r0 <- rv(); p0 <- r0$prod; e0 <- r0$ej
  if ((p0 == "vanilla" && e0 == "american") || (p0 == "barrier" && e0 == "american")) {
    return(div(class = "warn", strong("Monte Carlo not displayed for American options."), br(), "Basic Monte Carlo simulates paths forward and estimates discounted expected payoffs, but it does not directly solve the optimal early-exercise decision. The binomial tree is used because it compares continuation value and exercise value at each node."))
  }
  
  req(r0$mc)
  ic <- r0$ic
  ref <- r0$mc_ref
  
  mc_settings <- div(class = "rbox", h4("Monte Carlo parameters"),
                     fluidRow(
                       column(3, numericInput("mc_n_sims", "Number of simulations", value = isolate(r0$mc_n_sims %||% 10000), min = 100, max = 200000, step = 1000)),
                       column(3, numericInput("mc_n_steps", "Steps per path", value = isolate(r0$mc_n_steps %||% 252), min = 10, max = 2520, step = 10)),
                       column(3, numericInput("mc_semilla", "Random seed (0 = off)", value = isolate(r0$mc_semilla %||% 111), min = 0, step = 1)),
                       column(3, br(), actionButton("update_mc", "Update Monte Carlo", class = "btn-primary", style = "width:100%; margin-top:5px;"))),
                     p(class = "ic",
                       "Monte Carlo parameters are applied only after pressing Update Monte Carlo.", br(), "For vanilla options, pricing uses one exact jump to maturity; the steps input is used for visual paths and Euler comparison. For path-dependent options, steps per path determines how many times each simulated path is observed."))
  
  nota <- if (p0 == "asian") {
    div(class = "info", "Asian options are path-dependent because the payoff depends on the average price along the path.")
  } else if (p0 == "lookback") {
    div(class = "info", "Lookback options are path-dependent because the payoff depends on the maximum or minimum reached along the path.")
  } else if (p0 == "barrier") {
    div(class = "info", "Barrier options are path-dependent because the payoff depends on whether the barrier has been touched.")
  } else {
    NULL
  }
  
  stats_df <- data.frame(
    Statistic = c("Mean payoff", "Std dev payoff", "% ITM", "Max payoff", "95% CI width"),
    Value = round(
      c(
        mean(r0$mc$payoffs),
        sd(r0$mc$payoffs),
        100 * mean(r0$mc$payoffs > 0),
        max(r0$mc$payoffs),
        ic$ic_sup - ic$ic_inf),4))
  
  table_row <- if (p0 == "vanilla" && !is.null(r0$euler_df)) {
    fluidRow(
      column(6, h4("Payoff statistics"), div(class = "rbox", style = "height:240px;", clean_table(stats_df))),
      column(6, h4("GBM exact vs Euler-Maruyama"), div(class = "rbox", style = "height:240px;", clean_table(r0$euler_df))))
  } else {
    fluidRow(column(5, h4("Payoff statistics"), div(class = "rbox", clean_table(stats_df))))
  }
  
  tagList(mc_settings,nota,
          fluidRow(
            column(4,div(
              class = "rbox metric",
              p(class = "plabel", paste0("Monte Carlo price (M = ", format(ic$M, big.mark = ","), ")")),
              p(class = "pbig", format(round(ic$precio, 4), nsmall = 4)),
              p(class = "ic", paste0("95% CI: [", round(ic$ic_inf, 4), ", ", round(ic$ic_sup, 4), "]")))),
            column(4,div(class = "rbox metric",
                         p(class = "plabel", "Analytical reference"),
                         if (!is.null(ref$precio) && !is.na(ref$precio)) {
                           tagList(
                             p(class = "pbig", format(round(ref$precio, 4), nsmall = 4)),
                             p(class = "ic", ref$etiqueta)
                           )
                         } else {
                           p(style = "color:#697586", "Not available")
                         })),
            column(4,
                   div(class = "rbox metric", p(class = "plabel", "Difference"),
                       if (!is.null(ref$precio) && !is.na(ref$precio)) {
                         dif <- ic$precio - ref$precio
                         p(class = "pbig", style = paste0("color:", if (abs(dif) < 0.10) "#16a34a" else "#f97316"),
                           format(round(dif, 4), nsmall = 4))
                       } else {
                         p(style = "color:#697586", "N/A")
                       }))),
          
          h4("Simulated paths"),
          div(class = "rbox", plotlyOutput("mc_caminos_plot", height = "400px")),
          table_row)
  })
  
  output$mc_caminos_plot <- renderPlotly({
    req(rv(), rv()$mc)
    r0 <- rv(); cam <- r0$mc$caminos
    if (is.null(cam) || ncol(cam) <= 2) {
      sem <- if ((r0$mc_semilla %||% 111) > 0) r0$mc_semilla else 111
      S_plot <- if (!is.null(r0$S_star)) r0$S_star else r0$S
      q_plot <- if (!is.null(r0$S_star)) 0 else r0$q
      res_v <- mc_simular_gbm(S_plot,r0$r,q_plot,r0$sigma,r0$T,100,min(r0$ic$M,100),sem)
      cam <- res_v$caminos
    }
    np <- min(nrow(cam),100); nst <- ncol(cam)-1; ts <- seq(0,r0$T,length.out=nst+1)
    df <- data.frame()
    for(i in 1:np) df <- rbind(df, data.frame(t=ts, S=cam[i,], path=as.character(i)))
    p <- plot_ly(df, x=~t, y=~S, split=~path, type="scatter", mode="lines",
                 line=list(color="rgba(37,99,235,0.35)", width=1), hoverinfo="none", showlegend=FALSE) %>%
      layout(xaxis=list(title="Time"), yaxis=list(title="Underlying price"), margin=list(l=50,r=20,b=50,t=20),
             shapes=list(list(type="line", x0=0, x1=r0$T, y0=r0$K, y1=r0$K, line=list(color="#ef4444", dash="dash"))))
    if (r0$prod=="barrier" && !is.null(r0$B)) {
      p <- layout(p, shapes=list(
        list(type="line", x0=0, x1=r0$T, y0=r0$K, y1=r0$K, line=list(color="#ef4444", dash="dash")),
        list(type="line", x0=0, x1=r0$T, y0=r0$B, y1=r0$B, line=list(color="#f97316", dash="dot"))
      ))
    }
    p %>% config(displaylogo=FALSE)
  })
  
  output$tab_mc_greeks <- renderUI({
    req(rv(), rv()$gr_mc)
    r0 <- rv(); g <- r0$gr_mc
    gc_box <- function(lbl, val, form) div(class="rbox", h5(lbl), p(style="font-size:20px;font-weight:bold;", if(is.na(val)) "N/A" else round(as.numeric(val),6)), p(style="font-size:11px;color:#697586;", form))
    df <- NULL
    comparison <- NULL
    if (r0$prod=="vanilla" && r0$ej=="european") {
      gr_bs <- r0$gr_bs
      df <- data.frame(Greek=c("Delta","Gamma","Theta/day","Vega/1% sigma","Rho/1% r"), MC=round(c(g$delta,g$gamma,g$theta_diario,g$vega_pct,g$rho_pct),6), Black_Scholes=round(c(gr_bs$delta,gr_bs$gamma,gr_bs$theta_diario,gr_bs$vega_pct,gr_bs$rho_pct),6))
      df$Difference <- round(df$MC - df$Black_Scholes, 6)
      comparison <- tagList(h4("Comparison with Black-Scholes"), div(class="rbox", clean_table(df)))
    }
    tagList(
      h4("Monte Carlo Greeks"),
      div(class="info", "Monte Carlo Greeks use up to 3,000 simulations per perturbation for computational efficiency. If the selected number of simulations is below 3,000, that value is used instead. The same random seed is reused across perturbations to reduce simulation noise."),
      fluidRow(
        column(3, gc_box("Delta", g$delta, "(V(S+h)-V(S-h))/(2h), h=1% S")),
        column(3, gc_box("Gamma", g$gamma, "(V(S+h)-2V+V(S-h))/h^2")),
        column(3, gc_box("Theta/day", g$theta_diario, "V(T-1/365)-V(T)")),
        column(3, gc_box("Vega", g$vega_pct, "Sensitivity to a 1% volatility change"))
      ),
      fluidRow(column(3, gc_box("Rho", g$rho_pct, "Sensitivity to a 1% rate change"))),
      comparison
    )
  })
}


############################################################################### 
# 5. LANZAMIENTO DE LA APLICACIÓN
############################################################################### 

shinyApp(ui = ui, server = server)