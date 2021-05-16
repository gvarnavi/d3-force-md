import {quadtree} from "d3-quadtree";
import constant from "./constant.js";
import jiggle from "./jiggle.js";
import {x, y} from "./simulation.js";

export default function() {
  var nodes,
      node,
      random,
      alpha,
      strength = constant(1), // strength now serves as a counter
      strengths,
      distanceMin2 = 0.75,
      distanceMax2 = Infinity,
      theta2 = 0.0625, // lennard-jones type interactions mostly short-range
      N=12, M=6; // lennard-jones exponents

  function force(_) {
    var i, n = nodes.length, tree = quadtree(nodes, x, y).visitAfter(accumulate);
    for (alpha = _, i = 0; i < n; ++i) {
      node = nodes[i];
      node.energy = -1; //self-energy
      node.force_x = 0, node.force_y =0;
      tree.visit(apply);
    }
  }

  function initialize() {
    if (!nodes) return;
    var i, n = nodes.length, node;
    strengths = new Array(n);
    for (i = 0; i < n; ++i) node = nodes[i], strengths[node.index] = +strength(node, i, nodes);
    if (M==N) N += 1; //avoid dividing by zero
  }

  function accumulate(quad) {
    var strength = 0, q, c, weight = 0, x, y, i;

    // For internal nodes, accumulate forces from child quadrants.
    if (quad.length) {
      for (x = y = i = 0; i < 4; ++i) {
        if ((q = quad[i]) && (c = Math.abs(q.value))) {
          strength += q.value, weight += c, x += c * q.x, y += c * q.y;
        }
      }
      quad.x = x / weight;
      quad.y = y / weight;
    }

    // For leaf nodes, accumulate forces from coincident quadrants.
    else {
      q = quad;
      q.x = q.data.x;
      q.y = q.data.y;
      do strength += strengths[q.data.index];
      while (q = q.next);
    }

    quad.value = strength;
  }

  function apply(quad, x1, _, x2) {
    if (!quad.value) return true;

    var x = quad.x - node.x,
        y = quad.y - node.y,
        w = x2 - x1,
        l = x * x + y * y,
        force_prefactor;

    // Apply the Barnes-Hut approximation if possible.
    // Limit forces for very close nodes; randomize direction if coincident.
    if (w * w / theta2 < l) {
      if (l < distanceMax2) {
        if (x === 0) x = jiggle(random), l += x * x;
        if (y === 0) y = jiggle(random), l += y * y;
        if (l < distanceMin2) l = Math.sqrt(distanceMin2 * l);

        force_prefactor = M*N/(M-N)*(1-Math.pow(l,(N-M)/2))/Math.pow(l,(N+2)/2);

        node.energy += (N*Math.pow(l,-M/2)-M*Math.pow(l,-N/2))/(M-N)*quad.value;
        node.force_x += force_prefactor*x*quad.value*alpha;
        node.force_y += force_prefactor*y*quad.value*alpha;

      }
      return true;
    }

    // Otherwise, process points directly.
    else if (quad.length || l >= distanceMax2) return;

    // Limit forces for very close nodes; randomize direction if coincident.
    if (quad.data !== node || quad.next) {
      if (x === 0) x = jiggle(random), l += x * x;
      if (y === 0) y = jiggle(random), l += y * y;
      if (l < distanceMin2) l = Math.sqrt(distanceMin2 * l);
    }

    force_prefactor = M*N/(M-N)*(1-Math.pow(l,(N-M)/2))/Math.pow(l,(N+2)/2);

    do if (quad.data !== node) {
      w = strengths[quad.data.index];
      node.energy += (N*Math.pow(l,-M/2)-M*Math.pow(l,-N/2))/(M-N)*w;
      node.force_x += force_prefactor*x*w*alpha;
      node.force_y += force_prefactor*y*w*alpha;
    } while (quad = quad.next);
  }

  force.initialize = function(_nodes, _random) {
    nodes = _nodes;
    random = _random;
    initialize();
  };

  force.strength = function(_) {
    return arguments.length ? (strength = typeof _ === "function" ? _ : constant(+_), initialize(), force) : strength;
  };

  force.distanceMin = function(_) {
    return arguments.length ? (distanceMin2 = _ * _, force) : Math.sqrt(distanceMin2);
  };

  force.distanceMax = function(_) {
    return arguments.length ? (distanceMax2 = _ * _, force) : Math.sqrt(distanceMax2);
  };

  force.theta = function(_) {
    return arguments.length ? (theta2 = _ * _, force) : Math.sqrt(theta2);
  };
  
  force.repulsivePower = function(_) {
    return arguments.length ? (N = _, force) : N;
  };

  force.attractivePower = function(_) {
    return arguments.length ? (M = _, force) : M;
  };

  return force;
}
