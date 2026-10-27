/*
	Dopetrope by HTML5 UP
	html5up.net | @ajlkn
	Free for personal and commercial use under the CCA 3.0 license (html5up.net/license)
*/

(function($) {

	var	$window = $(window),
		$body = $('body');

	// Breakpoints.
		breakpoints({
			xlarge:  [ '1281px',  '1680px' ],
			large:   [ '981px',   '1280px' ],
			medium:  [ '737px',   '980px'  ],
			small:   [ null,      '736px'  ]
		});

	// Play initial animations on page load.
		$window.on('load', function() {
			window.setTimeout(function() {
				$body.removeClass('is-preload');
			}, 100);
		});

	// Dropdowns.
		$('#nav > ul').dropotron({
			mode: 'fade',
			noOpenerFade: true,
			alignment: 'center'
		});

	// Nav.

		// Title Bar.
			$(
				'<div id="titleBar">' +
					'<a href="#navPanel" class="toggle"></a>' +
				'</div>'
			)
				.appendTo($body);

		// Panel.
			$(
				'<div id="navPanel">' +
					'<nav>' +
						$('#nav').navList() +
					'</nav>' +
				'</div>'
			)
				.appendTo($body)
				.panel({
					delay: 500,
					hideOnClick: true,
					hideOnSwipe: true,
					resetScroll: true,
					resetForms: true,
					side: 'left',
					target: $body,
					visibleClass: 'navPanel-visible'
				});

	// Video: pause button, and don't autoplay for people who prefer reduced motion.
		$('video.intro-video').each(function() {
			var v = this, $btn = $(v).siblings('.video-toggle');
			function sync() {
				$btn.text(v.paused ? 'Play animation' : 'Pause animation')
					.attr('aria-pressed', v.paused ? 'true' : 'false');
			}
			if (window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
				v.removeAttribute('autoplay');
				v.pause();
			}
			$btn.on('click', function() { if (v.paused) { v.play(); } else { v.pause(); } });
			$(v).on('play pause', sync);
			sync();
		});

})(jQuery);
