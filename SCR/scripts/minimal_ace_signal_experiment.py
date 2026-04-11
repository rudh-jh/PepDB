# -*- coding: utf-8 -*-
"""
PepDB 最小可运行实验：
用 ACE 高置信二三肽稳定序列表，做活性分层分类
特征对比：AAC vs Signal vs Fusion
模型：LinearSVC
评估：Stratified K-Fold + Accuracy / Macro-P / Macro-R / Macro-F1 / 混淆矩阵
"""

from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")   # 强制使用非交互式后端，适合保存图片
import matplotlib.pyplot as plt


plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.svm import LinearSVC
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import (
    accuracy_score,
    precision_recall_fscore_support,
    confusion_matrix,
    ConfusionMatrixDisplay,
)

# =========================
# 1. 路径配置
# =========================

REPO_ROOT = Path(__file__).resolve().parents[2]

DATA_FILE = REPO_ROOT / "DB" / "worksets" / "ace" / "ahtpdb" / "高置信短肽工作集" / "高置信二三肽稳定序列表.csv"
OUT_DIR = REPO_ROOT / "DB" / "analysis" / "ace" / "minimal_signal_experiment"

# =========================
# 2. 氨基酸理化属性表
#    尽量轻量，够课程实验用
# =========================

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

HYDRO = {  # Kyte-Doolittle hydrophobicity
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}

CHARGE = {  # 简化净电荷
    "A": 0.0, "C": 0.0, "D": -1.0, "E": -1.0, "F": 0.0,
    "G": 0.0, "H": 0.1, "I": 0.0, "K": 1.0, "L": 0.0,
    "M": 0.0, "N": 0.0, "P": 0.0, "Q": 0.0, "R": 1.0,
    "S": 0.0, "T": 0.0, "V": 0.0, "W": 0.0, "Y": 0.0,
}

POLARITY = {
    "A": 8.1, "C": 5.5, "D": 13.0, "E": 12.3, "F": 5.2,
    "G": 9.0, "H": 10.4, "I": 5.2, "K": 11.3, "L": 4.9,
    "M": 5.7, "N": 11.6, "P": 8.0, "Q": 10.5, "R": 10.5,
    "S": 9.2, "T": 8.6, "V": 5.9, "W": 5.4, "Y": 6.2,
}

MW = {
    "A": 89.09, "C": 121.15, "D": 133.10, "E": 147.13, "F": 165.19,
    "G": 75.07, "H": 155.16, "I": 131.17, "K": 146.19, "L": 131.17,
    "M": 149.21, "N": 132.12, "P": 115.13, "Q": 146.15, "R": 174.20,
    "S": 105.09, "T": 119.12, "V": 117.15, "W": 204.23, "Y": 181.19,
}

PROPERTY_MAPS: Dict[str, Dict[str, float]] = {
    "hydro": HYDRO,
    "charge": CHARGE,
    "polarity": POLARITY,
    "mw": MW,
}

# =========================
# 3. 数据读取
# =========================

def robust_read_csv(path: Path) -> pd.DataFrame:
    """
    尽量兼容不同编码/换行格式。
    """
    if not path.exists():
        raise FileNotFoundError(f"找不到数据文件：{path}")

    # 尝试直接读取
    for enc in ["utf-8-sig", "utf-8", "gbk"]:
        try:
            df = pd.read_csv(path, encoding=enc)
            if df.shape[1] > 1:
                return df
        except Exception:
            pass

    # 如果直接读取失败，再做文本层面修复
    text = None
    for enc in ["utf-8-sig", "utf-8", "gbk"]:
        try:
            text = path.read_text(encoding=enc)
            break
        except Exception:
            continue

    if text is None:
        raise RuntimeError("无法读取 CSV 文件，请检查编码。")

    # 兼容奇怪换行
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    df = pd.read_csv(io.StringIO(text))
    return df


# =========================
# 4. 特征工程
# =========================

def clean_sequence(seq: str) -> str:
    seq = str(seq).strip().upper()
    seq = "".join([aa for aa in seq if aa in AA_LIST])
    return seq

def aac_features(seq: str) -> np.ndarray:
    seq = clean_sequence(seq)
    n = len(seq)
    if n == 0:
        return np.zeros(len(AA_LIST) + 1, dtype=float)

    feats = [seq.count(aa) / n for aa in AA_LIST]
    feats.append(float(n))  # 肽长
    return np.array(feats, dtype=float)

def safe_autocorr(x: np.ndarray, lag: int) -> float:
    if len(x) <= lag or np.std(x) < 1e-12:
        return 0.0
    x1 = x[:-lag]
    x2 = x[lag:]
    if np.std(x1) < 1e-12 or np.std(x2) < 1e-12:
        return 0.0
    return float(np.corrcoef(x1, x2)[0, 1])

def signal_features_for_property(seq: str, prop_map: Dict[str, float]) -> List[float]:
    seq = clean_sequence(seq)
    x = np.array([prop_map[aa] for aa in seq], dtype=float)
    if len(x) == 0:
        return [0.0] * 10

    # 时域
    mean_ = float(np.mean(x))
    std_ = float(np.std(x))
    min_ = float(np.min(x))
    max_ = float(np.max(x))
    energy_ = float(np.mean(x ** 2))

    # 相关域
    ac1 = safe_autocorr(x, 1)
    ac2 = safe_autocorr(x, 2)

    # 频域
    fft_vals = np.fft.rfft(x)
    mag = np.abs(fft_vals)
    fft_energy = float(np.sum(mag ** 2))

    if np.sum(mag) < 1e-12:
        centroid = 0.0
    else:
        centroid = float(np.sum(np.arange(len(mag)) * mag) / np.sum(mag))

    dom_freq = float(np.argmax(mag)) if len(mag) > 0 else 0.0
    dom_freq_norm = dom_freq / max(1, len(mag) - 1)

    length_ = float(len(x))

    return [
        mean_, std_, min_, max_, energy_,
        ac1, ac2,
        fft_energy, centroid, dom_freq_norm + length_ * 0.0  # 保持10维，避免花哨
    ]

def signal_features(seq: str) -> np.ndarray:
    feats = []
    for _, prop_map in PROPERTY_MAPS.items():
        feats.extend(signal_features_for_property(seq, prop_map))
    # 加一个长度特征
    feats.append(float(len(clean_sequence(seq))))
    return np.array(feats, dtype=float)

def build_feature_matrix(seqs: List[str], mode: str) -> np.ndarray:
    if mode == "aac":
        return np.vstack([aac_features(s) for s in seqs])
    elif mode == "signal":
        return np.vstack([signal_features(s) for s in seqs])
    elif mode == "fusion":
        return np.vstack([
            np.concatenate([aac_features(s), signal_features(s)]) for s in seqs
        ])
    else:
        raise ValueError(f"未知 mode: {mode}")


# =========================
# 5. 评估
# =========================

def evaluate_feature_set(
    X: np.ndarray,
    y: np.ndarray,
    class_names: List[str],
    out_prefix: Path,
) -> Dict[str, float]:
    class_counts = pd.Series(y).value_counts()
    min_count = int(class_counts.min())
    n_splits = max(2, min(5, min_count))

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    model = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LinearSVC(class_weight="balanced", random_state=42, dual="auto")),
    ])

    y_pred = cross_val_predict(model, X, y, cv=cv)

    acc = accuracy_score(y, y_pred)
    p, r, f1, _ = precision_recall_fscore_support(
        y, y_pred, average="macro", zero_division=0
    )

    cm = confusion_matrix(y, y_pred)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=class_names)

    fig, ax = plt.subplots(figsize=(6, 5))
    disp.plot(ax=ax, cmap="Blues", colorbar=False)
    ax.set_title(f"{out_prefix.stem} 混淆矩阵")
    fig.tight_layout()
    fig.savefig(out_prefix.with_suffix(".png"), dpi=200)
    plt.close(fig)

    return {
        "accuracy": round(float(acc), 4),
        "macro_precision": round(float(p), 4),
        "macro_recall": round(float(r), 4),
        "macro_f1": round(float(f1), 4),
        "n_splits": n_splits,
    }


# =========================
# 6. 主流程
# =========================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = robust_read_csv(DATA_FILE)

    # 只保留最关键列
    seq_col = "序列"
    label_col = "活性分层建议"

    if seq_col not in df.columns or label_col not in df.columns:
        raise KeyError(
            f"数据里没找到关键列。当前列名：{list(df.columns)}"
        )

    df = df[[seq_col, label_col]].copy()
    df[seq_col] = df[seq_col].astype(str).map(clean_sequence)
    df[label_col] = df[label_col].astype(str).str.strip()

    # 去掉空序列和空标签
    df = df[(df[seq_col] != "") & (df[label_col] != "")].drop_duplicates()

    # 去掉样本数过少的类别，避免 CV 崩
    vc = df[label_col].value_counts()
    keep_labels = vc[vc >= 2].index
    df = df[df[label_col].isin(keep_labels)].copy()

    print("===== 数据概况 =====")
    print(f"文件: {DATA_FILE}")
    print(f"样本数: {len(df)}")
    print("类别分布:")
    print(df[label_col].value_counts())
    print()

    seqs = df[seq_col].tolist()
    labels = df[label_col].tolist()

    le = LabelEncoder()
    y = le.fit_transform(labels)
    class_names = list(le.classes_)

    results = []
    for mode in ["aac", "signal", "fusion"]:
        print(f"===== 开始评估: {mode} =====")
        X = build_feature_matrix(seqs, mode)
        metrics = evaluate_feature_set(
            X=X,
            y=y,
            class_names=class_names,
            out_prefix=OUT_DIR / f"confusion_{mode}",
        )
        metrics["feature_set"] = mode
        results.append(metrics)
        print(metrics)
        print()

    result_df = pd.DataFrame(results)[
        ["feature_set", "accuracy", "macro_precision", "macro_recall", "macro_f1", "n_splits"]
    ]
    result_df.to_csv(OUT_DIR / "metrics_summary.csv", index=False, encoding="utf-8-sig")

    with open(OUT_DIR / "run_summary.json", "w", encoding="utf-8") as f:
        json.dump(
            {
                "data_file": str(DATA_FILE),
                "n_samples": int(len(df)),
                "class_distribution": df[label_col].value_counts().to_dict(),
                "results": results,
            },
            f,
            ensure_ascii=False,
            indent=2,
        )

    print("===== 全部完成 =====")
    print(result_df)
    print(f"\n结果已保存到: {OUT_DIR}")


if __name__ == "__main__":
    main()