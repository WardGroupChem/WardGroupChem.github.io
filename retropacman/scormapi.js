/* Minimal SCORM 1.2 API wrapper.
   Finds the LMS API object in a parent/opener window (as Moodle's SCORM
   player provides) and exposes simple init/setScore/finish helpers.
   If no LMS API is found (e.g. testing the game outside Moodle), it
   silently no-ops so the game still works standalone. */
(function (global) {
  var API = null;
  var initialized = false;

  function findAPI(win) {
    var attempts = 0;
    while (win.API == null && win.parent != null && win.parent !== win && attempts < 10) {
      attempts++;
      win = win.parent;
    }
    return win.API || null;
  }

  function locateAPI() {
    var theAPI = null;
    if (window.API) {
      theAPI = window.API;
    } else if (window.parent && window.parent !== window) {
      theAPI = findAPI(window.parent);
    }
    if (!theAPI && window.opener) {
      theAPI = findAPI(window.opener);
    }
    return theAPI;
  }

  var SCORM = {
    ready: false,
    init: function () {
      API = locateAPI();
      if (!API) {
        console.warn('SCORM API not found — running in standalone/test mode.');
        return false;
      }
      var result = API.LMSInitialize('');
      initialized = (result === 'true' || result === true);
      this.ready = initialized;
      if (initialized) {
        // Resume in-progress status if not already completed
        try {
          var status = API.LMSGetValue('cmi.core.lesson_status');
          if (status === 'not attempted' || status === '') {
            API.LMSSetValue('cmi.core.lesson_status', 'incomplete');
            API.LMSCommit('');
          }
        } catch (e) {}
      }
      return initialized;
    },
    setScore: function (rawScore, maxScore) {
      if (!API || !initialized) return;
      try {
        API.LMSSetValue('cmi.core.score.raw', String(rawScore));
        API.LMSSetValue('cmi.core.score.min', '0');
        API.LMSSetValue('cmi.core.score.max', String(maxScore));
        API.LMSCommit('');
      } catch (e) { console.warn('SCORM setScore failed', e); }
    },
    setComplete: function (passed) {
      if (!API || !initialized) return;
      try {
        API.LMSSetValue('cmi.core.lesson_status', passed ? 'passed' : 'failed');
        API.LMSCommit('');
      } catch (e) { console.warn('SCORM setComplete failed', e); }
    },
    finish: function () {
      if (!API || !initialized) return;
      try {
        API.LMSCommit('');
        API.LMSFinish('');
      } catch (e) { console.warn('SCORM finish failed', e); }
    }
  };

  global.SCORM = SCORM;
})(window);
