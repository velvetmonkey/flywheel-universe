# Flywheel Universe — Live Demo

This directory backs the GitHub Pages site at
<https://velvetmonkey.github.io/flywheel-universe/>.

`index.html` is a deploy copy of [`demos/hebbian-kuramoto.html`](../demos/hebbian-kuramoto.html) —
the canonical source lives there. Update the demo in `demos/`, then sync
the copy here when cutting a release:

```sh
cp demos/hebbian-kuramoto.html docs/index.html
git commit -am "docs: sync demo"
```

## What's the demo?

An interactive Hebbian-Kuramoto network — random clocks with adaptive
coupling that learn each other. Knock some out and watch the network
remember (or forget). Five topologies, speed control, senescence and
learning toggles, perturbation slider.
