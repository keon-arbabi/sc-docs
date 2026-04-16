// Animated benchmark bar chart — grows bars on scroll into view
(function() {
  const data = {
    'Basic Workflow': {
      brisc: 427, scanpy: 27518, seurat: 43301
    },
    'Label Transfer': {
      brisc: 144, scanpy: 4367, seurat: 53217
    },
    'Pseudobulk DE': {
      brisc: 81, scanpy: 15256, seurat: 4328
    }
  };

  function fmt(s) {
    if (s >= 3600) return (s / 3600).toFixed(1) + 'h';
    if (s >= 60) return (s / 60).toFixed(1) + 'm';
    return s.toFixed(0) + 's';
  }

  function build() {
    const container = document.getElementById('benchmark-chart');
    if (!container) return;

    Object.entries(data).forEach(([workflow, times]) => {
      const group = document.createElement('div');
      group.className = 'bench-group';

      const title = document.createElement('div');
      title.className = 'bench-title';
      title.textContent = workflow;
      group.appendChild(title);

      const maxTime = Math.max(...Object.values(times));

      const libs = [
        { name: 'brisc', time: times.brisc, cls: 'bar-brisc' },
        { name: 'Scanpy', time: times.scanpy, cls: 'bar-scanpy' },
        { name: 'Seurat + BPCells', time: times.seurat, cls: 'bar-seurat' }
      ];

      libs.forEach(lib => {
        const row = document.createElement('div');
        row.className = 'bench-row';

        const label = document.createElement('div');
        label.className = 'bench-label';
        label.textContent = lib.name;

        const track = document.createElement('div');
        track.className = 'bench-track';

        const bar = document.createElement('div');
        bar.className = 'bench-bar ' + lib.cls;
        bar.dataset.pct = ((lib.time / maxTime) * 100).toFixed(1);
        bar.style.width = '0%';

        const val = document.createElement('span');
        val.className = 'bench-val';
        val.textContent = fmt(lib.time);

        track.appendChild(bar);
        row.appendChild(label);
        row.appendChild(track);
        row.appendChild(val);
        group.appendChild(row);
      });

      container.appendChild(group);
    });

    // Animate on scroll into view
    const observer = new IntersectionObserver((entries) => {
      entries.forEach(entry => {
        if (entry.isIntersecting) {
          const bars = container.querySelectorAll('.bench-bar');
          bars.forEach((bar, i) => {
            setTimeout(() => {
              bar.style.width = bar.dataset.pct + '%';
            }, i * 80);
          });
          observer.disconnect();
        }
      });
    }, { threshold: 0.2 });

    observer.observe(container);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', build);
  } else {
    build();
  }
})();
