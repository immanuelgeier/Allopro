import os
import matplotlib.pyplot as plt
import pandas as pd

# =============================================================================
# KONFIGURATION
# =============================================================================
FILE_NAME = "coarse_marginals_real.txt"  # Oder "coarse_marginals_shuffle.txt"
TARGET_OBJECT = 1  # Welches Objekt (Kugel i) soll geplottet werden?
# =============================================================================

# Prüfen, ob die Datei existiert
if not os.path.exists(FILE_NAME):
    print(
        f"Fehler: Die Datei '{FILE_NAME}' wurde nicht gefunden. "
        "Hast du den Fortran-Lauf schon gestartet?"
    )
    exit(1)

# Daten einlesen (Pandas ignoriert die Header-Zeile mit '#' automatisch als Kommentar)
df = pd.read_csv(FILE_NAME, sep=r"\s+", comment="#", names=["Pixel", "Object", "X_Future", "Y_Past_Tgt", "Z_Past_Src"])

# Filtern nach dem gewünschten Objekt
df_obj = df[df["Object"] == TARGET_OBJECT]

if df_obj.empty:
    print(
        f"Fehler: Objekt {TARGET_OBJECT} wurde in der Datei nicht gefunden. "
        f"Verfügbare Objekte: {df['Object'].unique()}"
    )
    exit(1)

# Plot erstellen
plt.figure(figsize=(10, 6))

# Balken leicht versetzen (wichtig bei nur 12 groben HEALPix-Pixeln für die Übersicht)
width = 0.25
pixel_indices = df_obj["Pixel"]

plt.bar(
    pixel_indices - width,
    df_obj["X_Future"],
    width=width,
    label="X (Future / Zukunft)",
    alpha=0.8,
)
plt.bar(
    pixel_indices,
    df_obj["Y_Past_Tgt"],
    width=width,
    label="Y (Past Tgt / Verg. Ziel)",
    alpha=0.8,
)
plt.bar(
    pixel_indices + width,
    df_obj["Z_Past_Src"],
    width=width,
    label="Z (Past Src / Verg. Quelle)",
    alpha=0.8,
)

# Layout & Beschriftung
plt.title(
    f"1D-Marginalien (Orientierungsverteilung) - Objekt {TARGET_OBJECT}\nDatei: {FILE_NAME}",
    fontsize=14,
    fontweight="bold",
)
plt.xlabel("HEALPix Pixel Index", fontsize=12)
plt.ylabel("Counts (Häufigkeit)", fontsize=12)

# Da es der grobe Healpix-Raum ist, macht es Sinn, die x-Achse auf die Pixel zu fixieren
plt.xticks(pixel_indices)
plt.grid(axis="y", linestyle="--", alpha=0.5)
plt.legend(fontsize=11)

plt.tight_layout()

# Plot anzeigen
plt.show()