import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "NIRNAYA — Sovereign Optimization Solver | SIH26119",
  description:
    "High-performance mathematical optimization solver for LP, MILP, Convex QP & MIQP. Problem SIH26119 | Team 151198 — Visionaries for Change | MRPL",
  keywords: ["optimization", "solver", "LP", "MILP", "QP", "MIQP", "SIH26119", "MRPL", "NIRNAYA"],
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en" suppressHydrationWarning>
      <head>
        <link rel="icon" href="data:image/svg+xml,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'><text y='.9em' font-size='90'>⚡</text></svg>" />
      </head>
      <body className="min-h-screen bg-[#07070f]">{children}</body>
    </html>
  );
}
