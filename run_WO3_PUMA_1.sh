#!/bin/bash
#SBATCH --partition=serial
#SBATCH --account=2026-i-fs0751
#SBATCH --qos=curso-cpu-long
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name="WO3_PUMA"
#SBATCH --mail-user=jose.alvarezcastrillo@ucr.ac.cr
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/dev/null
#SBATCH --error=slurm-%x-%j.err

# =============================================================
# USO:
#   sbatch run_WO3_PUMA.sh S2 1        ← Run 1 de S2 barrido amplio
#   sbatch run_WO3_PUMA.sh S2 2        ← Run 2 de S2 refinamiento
#   sbatch run_WO3_PUMA.sh S2 3        ← Run 3 de S2 ajuste fino
#
# Esperar el mail de END antes de lanzar el siguiente run.
# El -sol.txt del run anterior debe estar en SRC para INIT=9. 
# =============================================================

SAMPLE_KEY="${1}"
RUN_NUM="${2}"

if [[ -z "$SAMPLE_KEY" || -z "$RUN_NUM" ]]; then
    echo "ERROR: uso: sbatch run_WO3_PUMA.sh <S2|S5|S8> <1|2|3>"
    exit 1
fi

# ----------------------------------------------------------
# Parámetros por muestra
# ----------------------------------------------------------
case "$SAMPLE_KEY" in
    S2)
        FNAME="WO3_S2_10sccm"
        # Run 1
        DMIN_1="0065"; DMAX_1="0150"; DSTEP_1=5
        INFMIN_1=0350; INFMAX_1=0550; INFSTEP_1=25
        # Run 2 — estrecho alrededor de 95 nm, inflexión 400 nm
        QUAD_2="2.89e-05"
        DMIN_2="0088"; DMAX_2="0102"; DSTEP_2=1
        INFMIN_2=0390; INFMAX_2=0410; INFSTEP_2=5
        # Run 3 — fijo en 95 nm, inflexión 400 nm
        QUAD_3="<QUAD_RUN2>"
        DFIX_3="0095"; INFFIX_3=0400
        ;;
    S5)
        FNAME="WO3_S5_7sccm"
        # Run 1
        DMIN_1="0100"; DMAX_1="0160"; DSTEP_1=5
        INFMIN_1=0350; INFMAX_1=0550; INFSTEP_1=25
        # Run 2 — estrecho alrededor de 130 nm, inflexión 350 nm
        QUAD_2="6.17e-05"
        DMIN_2="0123"; DMAX_2="0137"; DSTEP_2=1
        INFMIN_2=0340; INFMAX_2=0360; INFSTEP_2=5
        # Run 3 — fijo en 130 nm, inflexión 350 nm
        QUAD_3="<QUAD_RUN2>"
        DFIX_3="0130"; INFFIX_3=0350
        ;;
    S8)
        FNAME="WO3_S8_5sccm"
        # Run 1
        DMIN_1="0100"; DMAX_1="0180"; DSTEP_1=5
        INFMIN_1=0300; INFMAX_1=0550; INFSTEP_1=25
        # Run 2 — actualizar tras Run 1
        QUAD_2="<QUAD_RUN1>"
        DMIN_2="<DMIN>"; DMAX_2="<DMAX>"; DSTEP_2=1
        INFMIN_2="<INFMIN>"; INFMAX_2="<INFMAX>"; INFSTEP_2=5
        # Run 3 — actualizar tras Run 2
        QUAD_3="<QUAD_RUN2>"
        DFIX_3="<D_RUN2>"; INFFIX_3="<INF_RUN2>"
        ;;
    *)
        echo "ERROR: muestra desconocida '$SAMPLE_KEY'. Usar S2, S5 o S8."
        exit 1
        ;;
esac

# ----------------------------------------------------------
# Parámetros PUMA comunes
# ----------------------------------------------------------
NLAYERS=4; SLAYER=2; SUBSTRATE=60; DATATYPE=T; NOBS=100
LMIN=0350; LMAX=1500

# Rangos de n y k (solo usados en Run 1)
N0INI=1.90; N0FIN=2.50; N0STEP=0.10
NFINI=1.80; NFFIN=2.20; NFSTEP=0.10
K0INI=0.00; K0FIN=0.30; K0STEP=0.02

# ----------------------------------------------------------
SRC="/home/jose.alvarezcastrillo/labavanzado/Lab_Avanzado"
TMPDIR=$(mktemp -d)
echo "Muestra: $FNAME  |  Run: $RUN_NUM"
echo "Corriendo PUMA en: $TMPDIR"

cp $SRC/puma_mod           $TMPDIR/
cp $SRC/${FNAME}-dat.txt   $TMPDIR/
cp $SRC/${FNAME}-sol.txt   $TMPDIR/ 2>/dev/null

cd $TMPDIR
rm -f $SRC/${FNAME}-inf.txt

# ----------------------------------------------------------
# RUN 1 — Barrido amplio (INIT=0)
# n y k: rangos amplios definidos arriba
# ----------------------------------------------------------
if [[ "$RUN_NUM" == "1" ]]; then
    echo "=== RUN 1 — Barrido amplio ==="
    ./puma_mod $FNAME \
        $NLAYERS $SLAYER $SUBSTRATE $DATATYPE $NOBS \
        $LMIN $LMAX 3000 1e+100 0 \
        $DMIN_1 $DMAX_1 $DSTEP_1 \
        $INFMIN_1 $INFMAX_1 $INFSTEP_1 \
        $N0INI $N0FIN $N0STEP \
        $NFINI $NFFIN $NFSTEP \
        $K0INI $K0FIN $K0STEP

# ----------------------------------------------------------
# RUN 2 — Refinamiento (INIT=9)
# Requiere -sol.txt del Run 1 en SRC
# n y k: NO se pasan, PUMA los toma del -sol.txt
# ----------------------------------------------------------
elif [[ "$RUN_NUM" == "2" ]]; then
    echo "=== RUN 2 — Refinamiento ==="
    ./puma_mod $FNAME \
        $NLAYERS $SLAYER $SUBSTRATE $DATATYPE $NOBS \
        $LMIN $LMAX 5000 $QUAD_2 9 \
        $DMIN_2 $DMAX_2 $DSTEP_2 \
        $INFMIN_2 $INFMAX_2 $INFSTEP_2

# ----------------------------------------------------------
# RUN 3 — Ajuste fino (INIT=9, espesor e inflexión fijos)
# Requiere -sol.txt del Run 2 en SRC
# Actualizar QUAD_3, DFIX_3, INFFIX_3 con los valores del Run 2
# ----------------------------------------------------------
elif [[ "$RUN_NUM" == "3" ]]; then
    echo "=== RUN 3 — Ajuste fino ==="
    ./puma_mod $FNAME \
        $NLAYERS $SLAYER $SUBSTRATE $DATATYPE $NOBS \
        $LMIN $LMAX 50000 $QUAD_3 9 \
        $DFIX_3 $DFIX_3 1 \
        $INFFIX_3 $INFFIX_3 1

else
    echo "ERROR: run desconocido '$RUN_NUM'. Usar 1, 2 o 3."
    exit 1
fi

echo "=== PUMA finalizado ==="

# ----------------------------------------------------------
# Copiar resultados a SRC y limpiar
# ----------------------------------------------------------
cp ${FNAME}-sol.txt ${FNAME}-inf.txt $SRC/ 2>/dev/null
echo "Resultados copiados a $SRC"

rm -rf $TMPDIR
echo "Carpeta temporal eliminada."