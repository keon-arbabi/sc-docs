window.BENCHMARK_DATA = {
  "subtitle": "192 CPUs, 755 GB RAM",
  "groups": {
    "Basic workflow": {
      "hardware": "cpu",
      "bars": {
        "brisc": 196.19,
        "scanpy": 18407.88,
        "seurat": 44255.71
      }
    },
    "Label transfer": {
      "hardware": "cpu",
      "bars": {
        "brisc": 47.13,
        "scanpy": 2192.61,
        "seurat": 49830.04
      }
    },
    "Pseudobulk differential expression": {
      "hardware": "cpu",
      "bars": {
        "brisc": 28.86,
        "scanpy": 2049.56,
        "seurat": 4255.12
      }
    },
    "Basic workflow · CPU vs GPU": {
      "hardware": "gpu",
      "note": "96 CPUs, 4× H100 GPU, 752 GB RAM",
      "bars": {
        "brisc": 276.74,
        "rapids": 546.23
      }
    }
  }
};
