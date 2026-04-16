// Override pydata-sphinx-theme's 3-state theme switcher to toggle only
// between light and dark (skip "auto" mode) with a single click.
(function() {
  function applyMode(mode) {
    document.documentElement.dataset.mode = mode;
    document.documentElement.dataset.theme = mode;
    try {
      localStorage.setItem('mode', mode);
      localStorage.setItem('theme', mode);
    } catch (e) {}
    document.querySelectorAll('.dropdown-menu').forEach(el => {
      if (mode === 'dark') el.classList.add('dropdown-menu-dark');
      else el.classList.remove('dropdown-menu-dark');
    });
  }

  function toggleTheme(e) {
    if (e) {
      e.preventDefault();
      e.stopImmediatePropagation();
    }
    const current = document.documentElement.dataset.theme || 'light';
    applyMode(current === 'dark' ? 'light' : 'dark');
  }

  function attach() {
    // If the auto mode is currently active, flip to light/dark first
    if (document.documentElement.dataset.mode === 'auto') {
      const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
      applyMode(prefersDark ? 'dark' : 'light');
    }
    // Replace click handler on every theme switch button
    document.querySelectorAll('.theme-switch-button').forEach(btn => {
      // Clone node to strip existing listeners, then wire ours
      const fresh = btn.cloneNode(true);
      btn.parentNode.replaceChild(fresh, btn);
      fresh.addEventListener('click', toggleTheme);
    });
  }

  // Close the mobile sidebar drawer when any nav link inside it is clicked
  // (including a link to the current page, which would otherwise leave the
  // drawer open since no navigation occurs).
  function wireSidebarDismissal() {
    const modals = ['#pst-primary-sidebar-modal', '#pst-secondary-sidebar-modal'];
    modals.forEach(sel => {
      const modal = document.querySelector(sel);
      if (!modal) return;
      modal.addEventListener('click', (e) => {
        const link = e.target.closest('a');
        if (link && modal.contains(link) && modal.open) {
          modal.close();
        }
      });
    });
  }

  // Run after pydata-sphinx-theme's own listeners are attached
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', () => setTimeout(() => {
      attach();
      wireSidebarDismissal();
    }, 0));
  } else {
    setTimeout(() => {
      attach();
      wireSidebarDismissal();
    }, 0);
  }
})();
