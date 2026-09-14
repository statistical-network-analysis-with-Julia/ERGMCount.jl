/* Add ecosystem links and accessibility enhancements without replacing native
 * Documenter controls or events. The footer also works without JavaScript.
 */
(() => {
    "use strict";

    function addEcosystemNavigation() {
        // Keep absolute ecosystem links on this origin in the combined preview.
        // Hash-only Documenter permalinks retain their native click handling.
        document.querySelectorAll("a[href]").forEach((anchor) => {
            const href = anchor.getAttribute("href");
            if (
                !/^https?:\/\/statistical-network-analysis-with-julia\.github\.io(?:\/|$)/i.test(
                    href,
                )
            )
                return;
            const url = new URL(href);
            if (
                url.hostname.toLowerCase() ===
                "statistical-network-analysis-with-julia.github.io"
            ) {
                anchor.setAttribute(
                    "href",
                    url.pathname + url.search + url.hash,
                );
            }
        });

        const main = document.querySelector("#documenter .docs-main");
        const article = document.getElementById("documenter-page");
        if (!main || !article || main.querySelector(".snwj-ecosystem")) return;

        const skip = document.createElement("a");
        skip.className = "snwj-skip-link";
        skip.href = "#documenter-page";
        skip.textContent = "Skip to content";
        article.setAttribute("tabindex", "-1");
        document.body.prepend(skip);

        const nav = document.createElement("nav");
        nav.className = "snwj-ecosystem";
        nav.setAttribute("aria-label", "Ecosystem navigation");

        const brand = document.createElement("a");
        brand.className = "snwj-brand";
        brand.href = "/";
        const mark = document.createElement("span");
        mark.className = "snwj-brand-mark";
        mark.setAttribute("aria-hidden", "true");
        for (let i = 0; i < 3; i++) mark.append(document.createElement("i"));
        brand.append(
            mark,
            document.createTextNode("Statistical Network Analysis with Julia"),
        );
        nav.append(brand);

        const links = document.createElement("div");
        links.className = "snwj-ecosystem-links";
        for (const [label, href] of [
            ["Packages", "/packages/"],
            ["Get started", "/getting-started/"],
            ["Capabilities", "/capabilities/"],
        ]) {
            const link = document.createElement("a");
            link.href = href;
            link.textContent = label;
            links.append(link);
        }
        nav.append(links);
        main.prepend(nav);

        // Documenter renders the mobile-menu icon without an accessible label.
        // Leave its event handling intact while naming the existing control.
        const menu = document.getElementById("documenter-sidebar-button");
        if (menu && !menu.hasAttribute("aria-label")) {
            menu.setAttribute("aria-label", "Open documentation navigation");
        }
        const version = document.getElementById("documenter-version-selector");
        if (
            version &&
            !version.hasAttribute("aria-label") &&
            !version.hasAttribute("aria-labelledby")
        ) {
            version.setAttribute("aria-label", "Documentation version");
        }
    }

    function makeDocstringsAccessible() {
        for (const summary of document.querySelectorAll(
            "details.docstring > summary",
        )) {
            const links = Array.from(summary.querySelectorAll("a[href]"));
            if (links.length === 0) continue;
            const linkRow = document.createElement("div");
            linkRow.className = "snwj-docstring-links";
            for (const link of links) {
                const label = document.createElement("span");
                label.className = link.className;
                const name = link.textContent.trim();
                // Keep the binding's code markup in the disclosure label. Move
                // the original link (including its href and any id) outside the
                // interactive summary; do not duplicate any fragment targets.
                while (link.firstChild) label.append(link.firstChild);
                link.replaceWith(label);
                if (link.classList.contains("docstring-binding")) {
                    link.textContent = "Permalink";
                    link.setAttribute("aria-label", `Permalink to ${name}`);
                } else {
                    link.textContent = name;
                }
                linkRow.append(link);
            }
            summary.after(linkRow);
        }
    }

    function watchScrollableExamples() {
        const article = document.getElementById("documenter-page");
        if (!article || article.dataset.snwjScrollReady) return;
        article.dataset.snwjScrollReady = "true";
        const originalTabIndices = new WeakMap();
        let pending = false;

        function refresh() {
            pending = false;
            for (const node of article.querySelectorAll(
                "pre, pre > code, table, .math-container",
            )) {
                const style = getComputedStyle(node);
                const scrolls = (value) =>
                    value === "auto" || value === "scroll";
                const overflowing =
                    node.clientWidth > 0 &&
                    ((scrolls(style.overflowX) &&
                        node.scrollWidth > node.clientWidth + 1) ||
                        (scrolls(style.overflowY) &&
                            node.scrollHeight > node.clientHeight + 1));
                if (overflowing) {
                    if (!originalTabIndices.has(node)) {
                        originalTabIndices.set(
                            node,
                            node.getAttribute("tabindex"),
                        );
                    }
                    node.setAttribute("tabindex", "0");
                } else if (originalTabIndices.has(node)) {
                    const original = originalTabIndices.get(node);
                    if (original === null) node.removeAttribute("tabindex");
                    else node.setAttribute("tabindex", original);
                    originalTabIndices.delete(node);
                }
            }
        }

        function scheduleRefresh() {
            if (pending) return;
            pending = true;
            requestAnimationFrame(refresh);
        }

        // Highlighting and web-font loading can change scroll widths after DOM
        // readiness. Collapsed docstrings acquire dimensions only when opened.
        const observer = new MutationObserver(scheduleRefresh);
        observer.observe(article, {
            childList: true,
            subtree: true,
            attributes: true,
            attributeFilter: ["class", "style"],
        });
        observer.observe(document.documentElement, {
            attributes: true,
            attributeFilter: ["class"],
        });
        article.addEventListener("toggle", scheduleRefresh, true);
        window.addEventListener("resize", scheduleRefresh, { passive: true });
        window.addEventListener("load", scheduleRefresh, { once: true });
        if (document.fonts) document.fonts.ready.then(scheduleRefresh);
        scheduleRefresh();
    }

    function enhanceDocumentation() {
        addEcosystemNavigation();
        makeDocstringsAccessible();
        watchScrollableExamples();
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", enhanceDocumentation, {
            once: true,
        });
    } else {
        enhanceDocumentation();
    }
})();
