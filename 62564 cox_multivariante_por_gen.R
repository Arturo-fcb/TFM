# Bloque verificado por el autor: sustituye al modelo de Cox multivariante conjunto de una versión anterior del TFM por cuatro modelos independientes, uno por gen (cada uno ajustado por edad al diagnóstico y clasificación de alto riesgo). Se ejecuta a continuación del bucle de Cox univariante en 62564.R, sobre el mismo objeto surv_data.
# =============================================================================
# Regresion de Cox multivariante — un modelo independiente por gen
# GSE62564 (cohorte SEQC, n = 498)
#
# Este bloque SUSTITUYE al modelo conjunto original:
#   cox_multi <- coxph(Surv(os_days, os_event) ~ ARPP21 + DLG2 + AXL + MYCN +
#                       age_days + high_risk, data = surv_data)
# Va justo despues del bucle de Cox univariante en 62564.R, y usa el mismo
# objeto `surv_data` que ya se construye en ese script (no hace falta tocar
# nada de la carga de datos ni del preprocesado).
# =============================================================================
 
library(survival)
 
genes <- c("ARPP21", "DLG2", "AXL", "MYCN")
modelos <- list()
 
for (gen in genes) {
  cat("\n===== Modelo multivariante:", gen, "=====\n")
  form <- as.formula(paste0("Surv(os_days, os_event) ~ ", gen, " + age_days + high_risk"))
  m <- coxph(form, data = surv_data)
  print(summary(m))
  cat("Concordancia:", summary(m)$concordance[1], "\n")
  modelos[[gen]] <- m
 
  cat("Test de riesgos proporcionales (cox.zph):\n")
  print(cox.zph(m))
}
 
# ── Tabla resumen (HR, IC 95 %, p-valor y concordancia por gen) ─────────────
cat("\n\n===== RESUMEN: HR del gen en su propio modelo ajustado =====\n")
for (gen in genes) {
  m  <- modelos[[gen]]
  s  <- summary(m)
  co <- s$coefficients[gen, ]
  ci <- s$conf.int[gen, ]
  cat(sprintf(
    "%-8s HR = %.3f  IC95%% = [%.3f - %.3f]  p = %.4g  n = %d  eventos = %d  concordancia = %.3f\n",
    gen, co["exp(coef)"], ci["lower .95"], ci["upper .95"], co["Pr(>|z|)"],
    s$n, s$nevent, s$concordance[1]
  ))
}
 
# ── Nota de reproducibilidad ─────────────────────────────────────────────────
# Los resultados esperados (con este mismo surv_data, sesion limpia) son:
#   ARPP21: HR=0.882  IC95=[0.803-0.970]  p=0.009187  concordancia=0.842
#   DLG2:   HR=0.940  IC95=[0.833-1.061]  p=0.317     concordancia=0.826
#   AXL:    HR=0.870  IC95=[0.747-1.013]  p=0.07299   concordancia=0.838
#   MYCN:   HR=1.162  IC95=[1.063-1.270]  p=0.0009539 concordancia=0.850
# Si al ejecutar esto te salen numeros distintos, antes de tocar nada del TFM
# compara surv_data (dim, NAs, valores de los 4 genes) con el que se genera
# en el bloque de carga de datos de 62564.R — lo mas probable es un desajuste
# ahi, no en este bloque de Cox en si.
 
