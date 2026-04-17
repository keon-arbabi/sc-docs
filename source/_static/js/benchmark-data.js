window.BENCHMARK_DATA = {
  "subtitle": "192 CPUs, 755 GB RAM",
  "groups": {
    "Basic workflow": {
      "hardware": "cpu",
      "bars": {
        "brisc": 427.01,
        "scanpy": 27517.63,
        "seurat": 43300.98
      }
    },
    "Label transfer": {
      "hardware": "cpu",
      "bars": {
        "brisc": 143.63,
        "scanpy": 4366.99,
        "seurat": 53216.51
      }
    },
    "Pseudobulk differential expression": {
      "hardware": "cpu",
      "bars": {
        "brisc": 81.31,
        "scanpy": 3827.64,
        "seurat": 4328.49
      }
    },
    "Basic workflow · CPU vs GPU": {
      "hardware": "gpu",
      "note": "96 CPUs, 4× H100 GPU, 752 GB RAM",
      "bars": {
        "brisc": 491.67,
        "rapids": 1426.6
      }
    }
  }
};
