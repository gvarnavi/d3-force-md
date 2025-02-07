import {dispatch} from "d3-dispatch";
import {timer} from "d3-timer";
import lcg from "./lcg.js";

export function x(d) {
  return d.x;
}

export function y(d) {
  return d.y;
}

var initialRadius = 10,
    initialAngle = Math.PI * (3 - Math.sqrt(5));

export default function(nodes) {
  var simulation,
      dt = 1,		
      alpha = 1,
      alphaMin = 0.001,
      alphaDecay = 1 - Math.pow(alphaMin, 1 / 300),
      alphaTarget = 0,
      velocityDecay = 0.6,
      forces = new Map(),
      stepper = timer(step),
      event = dispatch("tick", "end"),
      random = lcg();

  if (nodes == null) nodes = [];

  function step() {
    tick(dt);
    event.call("tick", simulation);
    if (alpha < alphaMin) {
      stepper.stop();
      event.call("end", simulation);
    }
  }

  function tick(iterations) {
    var i, n = nodes.length, node;

    if (iterations === undefined) iterations = 1;

    for (var k = 0; k < iterations; ++k) {
      alpha += (alphaTarget - alpha) * alphaDecay;
      
      // update velocities (half-step)
      for (i = 0; i < n; ++i) {
        node = nodes[i];

        if (node.fx == null) node.vx += 0.5*node.force_x*dt/node.mass;
        else node.vx = 0, node.force_x = 0, node.x = node.fx;
        if (node.fy == null) node.vy += 0.5*node.force_y*dt/node.mass;
        else node.vy = 0, node.force_y = 0, node.y = node.fy;
      }
      
      // update positions
      for (i = 0; i < n; ++i) {
        node = nodes[i];
        
        if (node.fx == null) node.x += dt*(node.vx *= velocityDecay);
        else node.x = node.fx;
        if (node.fy == null) node.y += dt*(node.vy *= velocityDecay);
        else node.y = node.fy;
      }

      // compute new forces and energies based on positions
      forces.forEach(function(force) {
        force(alpha);
      });

      // update velocities (full-step)
      for (i = 0; i < n; ++i) {
        node = nodes[i];

        if (node.fx == null) node.vx += 0.5*node.force_x*dt/node.mass;
        else node.vx = 0, node.force_x = 0, node.x = node.fx;
        if (node.fy == null) node.vy += 0.5*node.force_y*dt/node.mass;
        else node.vy = 0, node.force_y = 0, node.y = node.fy;
      }

    }

    return simulation;
  }

  function initializeNodes() {
    for (var i = 0, n = nodes.length, node; i < n; ++i) {

      node = nodes[i], node.index = i;

      if (node.fx != null) node.x = node.fx;
      if (node.fy != null) node.y = node.fy;

      if (isNaN(node.x) || isNaN(node.y)) {
        var radius = initialRadius * Math.sqrt(0.5 + i), angle = i * initialAngle;
        node.x = radius * Math.cos(angle);
        node.y = radius * Math.sin(angle);
      }

      if (isNaN(node.vx) || isNaN(node.vy)) {
        node.vx = node.vy = 0;
      }

      if (isNaN(node.force_x) || isNaN(node.force_y)) {
        node.force_x = node.force_y = 0;
      }

      if (isNaN(node.mass)) node.mass = 1;
    }
  }

  function initializeForce(force) {
    if (force.initialize) force.initialize(nodes, random);
    return force;
  }

  initializeNodes();

  return simulation = {
    tick: tick,

    restart: function() {
      return stepper.restart(step), simulation;
    },

    stop: function() {
      return stepper.stop(), simulation;
    },

    nodes: function(_) {
      return arguments.length ? (nodes = _, initializeNodes(), forces.forEach(initializeForce), simulation) : nodes;
    },

    alpha: function(_) {
      return arguments.length ? (alpha = +_, simulation) : alpha;
    },

    alphaMin: function(_) {
      return arguments.length ? (alphaMin = +_, simulation) : alphaMin;
    },

    
    dt: function(_) {
      return arguments.length ? (dt = +_, simulation) : dt;
    },

    alphaDecay: function(_) {
      return arguments.length ? (alphaDecay = +_, simulation) : +alphaDecay;
    },

    alphaTarget: function(_) {
      return arguments.length ? (alphaTarget = +_, simulation) : alphaTarget;
    },

    velocityDecay: function(_) {
      return arguments.length ? (velocityDecay = 1 - _, simulation) : 1 - velocityDecay;
    },

    randomSource: function(_) {
      return arguments.length ? (random = _, forces.forEach(initializeForce), simulation) : random;
    },

    force: function(name, _) {
      return arguments.length > 1 ? ((_ == null ? forces.delete(name) : forces.set(name, initializeForce(_))), simulation) : forces.get(name);
    },

    find: function(x, y, radius) {
      var i = 0,
          n = nodes.length,
          dx,
          dy,
          d2,
          node,
          closest;

      if (radius == null) radius = Infinity;
      else radius *= radius;

      for (i = 0; i < n; ++i) {
        node = nodes[i];
        dx = x - node.x;
        dy = y - node.y;
        d2 = dx * dx + dy * dy;
        if (d2 < radius) closest = node, radius = d2;
      }

      return closest;
    },

    on: function(name, _) {
      return arguments.length > 1 ? (event.on(name, _), simulation) : event.on(name);
    }
  };
}
