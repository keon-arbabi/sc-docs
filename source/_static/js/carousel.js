// 3-card code carousel. One card is active (centered, full opacity), the
// other two flank it at reduced opacity. Arrows rotate; clicking a flanking
// card also activates it.
(function() {
  function apply(cards, activeIdx) {
    const n = cards.length;
    cards.forEach((c, i) => {
      c.classList.remove('is-active', 'is-prev', 'is-next');
      c.setAttribute('aria-hidden', i === activeIdx ? 'false' : 'true');
      if (i === activeIdx) c.classList.add('is-active');
      else if (i === (activeIdx - 1 + n) % n) c.classList.add('is-prev');
      else if (i === (activeIdx + 1) % n) c.classList.add('is-next');
    });
  }

  function rotate(carousel, dir) {
    const cards = [...carousel.querySelectorAll('.carousel-card')];
    const active = cards.findIndex(c => c.classList.contains('is-active'));
    const start = active < 0 ? 0 : active;
    const next = (start + dir + cards.length) % cards.length;
    apply(cards, next);
  }

  function activate(carousel, idx) {
    const cards = [...carousel.querySelectorAll('.carousel-card')];
    apply(cards, idx);
  }

  function init() {
    document.querySelectorAll('.code-carousel').forEach(carousel => {
      const leftBtn = carousel.querySelector('.carousel-arrow-left');
      const rightBtn = carousel.querySelector('.carousel-arrow-right');
      if (leftBtn) leftBtn.addEventListener('click', e => {
        e.preventDefault(); rotate(carousel, -1);
      });
      if (rightBtn) rightBtn.addEventListener('click', e => {
        e.preventDefault(); rotate(carousel, +1);
      });
      // Click on a flanking card promotes it to active.
      carousel.querySelectorAll('.carousel-card').forEach((card, i) => {
        card.addEventListener('click', e => {
          if (card.classList.contains('is-active')) return;
          // ignore clicks inside the active card's code
          e.preventDefault();
          activate(carousel, i);
        });
      });
      // Keyboard: left/right arrow keys when focused inside carousel.
      carousel.addEventListener('keydown', e => {
        if (e.key === 'ArrowLeft') { e.preventDefault(); rotate(carousel, -1); }
        else if (e.key === 'ArrowRight') { e.preventDefault(); rotate(carousel, +1); }
      });
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
