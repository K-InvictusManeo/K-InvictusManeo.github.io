document.documentElement.classList.add("js");

const scenes = Array.from(document.querySelectorAll("[data-video-scene]"));
const steps = Array.from(document.querySelectorAll("[data-video-step]"));
const playButton = document.querySelector("[data-video-play]");
const video = document.querySelector(".virtual-video");
const sceneDuration = 5800;
const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
let currentScene = 0;
let playbackTimer = null;

function showScene(index) {
  currentScene = (index + scenes.length) % scenes.length;
  scenes.forEach((scene, sceneIndex) => {
    const active = sceneIndex === currentScene;
    scene.classList.toggle("active", active);
    scene.hidden = !active;
  });
  steps.forEach((step, stepIndex) => {
    step.classList.toggle("active", stepIndex === currentScene);
    step.classList.toggle("complete", stepIndex < currentScene);
  });
}

function stopPlayback() {
  if (playbackTimer) window.clearInterval(playbackTimer);
  playbackTimer = null;
  video?.classList.remove("is-playing");
  playButton?.setAttribute("aria-pressed", "false");
  if (playButton) playButton.innerHTML = "<span>▶</span> 播放";
}

function startPlayback() {
  stopPlayback();
  playButton?.setAttribute("aria-pressed", "true");
  if (playButton) playButton.innerHTML = "<span>Ⅱ</span> 暂停";
  video?.classList.add("is-playing");
  playbackTimer = window.setInterval(() => showScene(currentScene + 1), sceneDuration);
}

steps.forEach((step) => {
  step.addEventListener("click", () => {
    stopPlayback();
    showScene(Number(step.dataset.videoStep));
  });
});

playButton?.addEventListener("click", () => {
  if (playbackTimer) stopPlayback();
  else startPlayback();
});

document.addEventListener("visibilitychange", () => {
  if (document.hidden) stopPlayback();
});

showScene(0);
if (reducedMotion) stopPlayback();
else startPlayback();

const revealSections = document.querySelectorAll(".reveal-section");
if (reducedMotion || !("IntersectionObserver" in window)) {
  revealSections.forEach((section) => section.classList.add("is-visible"));
} else {
  const revealObserver = new IntersectionObserver(
    (entries) => entries.forEach((entry) => {
      if (!entry.isIntersecting) return;
      entry.target.classList.add("is-visible");
      revealObserver.unobserve(entry.target);
    }),
    { threshold: 0.12 }
  );
  revealSections.forEach((section) => revealObserver.observe(section));
}
